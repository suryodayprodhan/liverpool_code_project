# Running the Polymer Intra-Chain Mobility Code on Windows

This code has two parts:

1. **Fortran modules** (`static_disordered_hamiltonian`, `calculate_ll`,
   `hopping_rate`, `steady_state_mobility`) that must be compiled with f2py
   into Python-importable modules.
2. **Python scripts** (`polymer_intra_chain_mobility.py`, `diffusivity.py`,
   `mobility_sweep.py`, etc.) that import those compiled modules.

The code is built with the **Intel oneAPI toolkit** (the `ifx` Fortran
compiler and Intel MKL for LAPACK/BLAS) for optimized, high-performance
numerical routines.

> **BUILD SYSTEM NOTE (Python 3.12+).** Python 3.12 removed the standard-
> library `distutils` module. As a result, f2py uses the **meson** build
> backend instead of the legacy distutils backend. The old `--fcompiler`
> flag is ignored by meson — compilers are selected via the `FC` and `CC`
> environment variables instead. The provided `build.bat` already handles
> this. You must install **meson** and **ninja** (via pip) before building.

There are two ways to run this on Windows. **Option A (WSL) is strongly
recommended** — it is far less painful than native Windows.

---

## OPTION A: With WSL (Recommended)

WSL (Windows Subsystem for Linux) runs a real Linux environment inside
Windows, where the Intel toolchain and f2py build behave exactly as on Linux.

### 1. Install WSL

Open PowerShell **as Administrator**:

    wsl --install

Restart when prompted (installs Ubuntu by default). Set a username and
password on first launch.

### 2. Follow the Linux instructions

Inside the Ubuntu terminal, follow the separate **Linux instructions** file.
In brief:

    # System tools
    sudo apt update
    sudo apt install -y build-essential python3 python3-pip python3-venv python3-dev

    # Intel oneAPI (Base + HPC toolkits)
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
        | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
        | sudo tee /etc/apt/sources.list.d/oneAPI.list
    sudo apt update
    sudo apt install -y intel-basekit intel-hpckit

    # Python environment
    python3 -m venv gnn
    source gnn/bin/activate
    pip install --upgrade pip
    pip install numpy==1.26.4 scipy matplotlib meson ninja

    # Build
    source /opt/intel/oneapi/setvars.sh
    chmod +x f2py_compile.sh
    ./f2py_compile.sh

    # Run
    python polymer_intra_chain_mobility.py

Your Windows files are visible under `/mnt/c/...` inside WSL, so you can keep
the project on your Windows drive.

---

## OPTION B: Native Windows (No WSL)

More involved, because Windows ships no Fortran compiler or LAPACK/BLAS. You
must install Intel oneAPI and use Python 3.12.

### 1. Install Python 3.12

Download **Python 3.12** (64-bit) from python.org. During install, check
**"Add python.exe to PATH"**.

### 2. Install Intel oneAPI for Windows

Download and install the **Intel oneAPI Base Toolkit** (provides MKL) and the
**Intel oneAPI HPC Toolkit** (provides the `ifx` Fortran compiler) from
Intel's website. After install, the environment is initialized with
`setvars.bat` (used in the build script below).

### 3. Set up the Python environment

Open Command Prompt (CMD):

    py -3.12 -m venv gnn
    gnn\Scripts\activate
    python --version

Confirm it says Python 3.12.x. Then install the required packages:

    pip install --upgrade pip
    pip install numpy scipy matplotlib meson ninja

Package notes:

- **numpy** — any recent version works (1.26.x or 2.x). The meson backend
  does not depend on `numpy.distutils` (which was removed in numpy 2.x),
  so numpy no longer needs to be pinned.
- **meson** and **ninja** — required by f2py's meson build backend on
  Python 3.12+. Without these, the build will fail.
- **scipy** is required by `diffusivity.py` (sparse matrix support).
- **matplotlib** is required by `mobility_sweep.py` for the plots.

### 4. Build the Fortran modules

Save the batch script (below) as `build.bat` in the project folder, then in
CMD (with the venv activated) run:

    build.bat

### 5. Run the code

    python polymer_intra_chain_mobility.py

---

## Running the Build Script on CMD

### Saving the file

1. Open a plain-text editor (Notepad, VS Code, or Notepad++).
2. Paste the batch-script contents (below).
3. Save as **`build.bat`** — in Notepad set "Save as type" to **"All Files"**
   (not "Text Documents"), or it becomes `build.bat.txt` and will not run.
4. Save it in the same folder as the `.f90` source files.

### Running the file

1. Open Command Prompt.
2. Navigate to the project folder:

       cd C:\Users\YourName\liverpool_code_project

3. Activate the environment:

       gnn\Scripts\activate

4. Run the script by typing its name (a `.bat` runs directly, no `python`):

       build.bat

On success you see "All modules built successfully." and new `.pyd` files
(the compiled modules) appear — the Windows equivalent of Linux `.so` files.

---

## The Windows Batch Build Script (build.bat)

On Python 3.12+, f2py uses the **meson** backend. Compilers are selected via
the `FC` and `CC` environment variables (not `--fcompiler`, which is ignored).
The script sets `FC=ifx` and `CC=icx` so meson picks up the Intel compilers
from oneAPI.

    @echo off
    setlocal enabledelayedexpansion

    REM ============================================================
    REM  Windows build script for the polymer mobility Fortran modules
    REM  Python 3.12+ with meson backend
    REM  Save as build.bat (NOT build.bat.txt)
    REM  Run from CMD with your Python venv activated:  build.bat
    REM ============================================================

    REM Adjust this path to match your oneAPI install location
    call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

    REM Tell meson which compilers to use (meson ignores --fcompiler)
    set FC=ifx
    set CC=icx

    set F2PY=python -m numpy.f2py -c
    set F90FLAGS=-O3 -xHost
    set LIBS=-lmkl_rt

    echo Building static_disordered_hamiltonian...
    %F2PY% %LIBS% datatype.f90 static_disordered_hamiltonian.f90 random_number_generator_fixedseed.f90 -m static_disordered_hamiltonian --f90flags="%F90FLAGS%" --verbose
    if errorlevel 1 goto error

    echo Building calculate_ll...
    %F2PY% %LIBS% datatype.f90 calculate_ll.f90 -m calculate_ll --f90flags="%F90FLAGS%" --verbose
    if errorlevel 1 goto error

    echo Building hopping_rate...
    %F2PY% %LIBS% datatype.f90 hopping_rate.f90 -m hopping_rate --f90flags="%F90FLAGS%" --verbose
    if errorlevel 1 goto error

    echo Building steady_state_mobility...
    %F2PY% %LIBS% datatype.f90 steady_state_mobility.f90 -m steady_state_mobility --f90flags="%F90FLAGS%" --verbose
    if errorlevel 1 goto error

    echo All modules built successfully.
    goto end

    :error
    echo Build failed.
    exit /b 1

    :end

Notes on the build script:

- **`set FC=ifx` / `set CC=icx`** — this is how meson discovers compilers on
  Python 3.12+. The old `--fcompiler=intelvem --f77exec=ifx --f90exec=ifx`
  flags are silently ignored by the meson backend and must not be relied on.
- **`-O3 -xHost`** — Intel optimization flags, passed via `--f90flags`. These
  are Intel-only; do not use them with gfortran. For debug builds, replace
  with `-check bounds -traceback -O0 -g`.
- **`-lmkl_rt`** — links Intel MKL for LAPACK/BLAS. Windows has no system
  `-llapack -lblas`, so MKL (bundled with oneAPI) is used instead.

---

## Troubleshooting

- **`.bat` opens in an editor instead of running** — it was saved as
  `build.bat.txt`. Rename to `build.bat`.
- **"'python' is not recognized"** — Python is not on PATH, or the venv is
  not activated. Activate with `gnn\Scripts\activate`, or use `py -3.12`.
- **`--fcompiler cannot be used with meson`** — you are using an old
  build script that passes `--fcompiler`. On Python 3.12+, remove that flag
  and set `FC=ifx` and `CC=icx` as environment variables instead (as the
  provided `build.bat` already does).
- **`Fortran compiler for the host machine: gfortran` in the meson output,
  even though you wanted Intel** — `FC` was not set before the build, so
  meson fell back to gfortran. Ensure `set FC=ifx` appears before the
  f2py calls, and that `setvars.bat` has been called so `ifx` is on PATH.
- **`gfortran: error: language Host not recognized`** — the Intel flag
  `-xHost` was passed to gfortran (see above — meson picked gfortran
  instead of ifx). Fix the `FC` variable.
- **`Cannot find -lmkl_rt`** — MKL not on the linker path. Ensure
  `setvars.bat` has been called (`echo %MKLROOT%` should show a path).
- **ifx not found** — `setvars.bat` was not called, or the path in the
  script is wrong. Check your oneAPI install location.
- **`No module named numpy`** — the venv is not activated, or numpy is not
  installed in it. Activate with `gnn\Scripts\activate` and reinstall.
- **`No module named mesonbuild`** — meson is not installed. Run
  `pip install meson ninja`.

When in doubt, use WSL (Option A). The native Windows Intel + f2py toolchain
is the main source of difficulty, and WSL avoids most of it.
