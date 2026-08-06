# Running the Polymer Intra-Chain Mobility Code on Windows

This code was developed on Linux. It has two parts:

1. **Fortran modules** (`static_disordered_hamiltonian`, `calculate_ll`,
   `hopping_rate`, `steady_state_mobility`) that must be compiled with f2py
   into Python-importable modules.
2. **Python scripts** (`polymer_intra_chain_mobility.py`, `diffusivity.py`,
   `mobility_sweep.py`, etc.) that import those compiled modules.

There are two ways to run it on Windows. **Option A (WSL) is strongly
recommended** — it is far less painful than native Windows.

---

## OPTION A: With WSL (Recommended)

WSL (Windows Subsystem for Linux) runs a real Linux environment inside
Windows. The original Linux build script works unchanged.

### 1. Install WSL

Open PowerShell **as Administrator** and run:

    wsl --install

Restart when prompted. This installs Ubuntu by default. Set a username and
password when it first launches.

### 2. Open the Ubuntu terminal

Search "Ubuntu" in the Start menu, or type `wsl` in a Command Prompt.

### 3. Install build tools, Python, and LAPACK/BLAS

    sudo apt update
    sudo apt install -y python3 python3-pip python3-venv gfortran liblapack-dev libblas-dev build-essential

### 4. (Optional) Install Intel oneAPI for ifx

Only needed if you specifically want the Intel Fortran compiler. Otherwise
skip this and use the gfortran build script (see step 6).

Follow Intel's instructions for installing oneAPI on Linux, then:

    source /opt/intel/oneapi/setvars.sh

### 5. Set up a Python environment

    python3 -m venv myenv
    source myenv/bin/activate
    pip install numpy==1.26.4 matplotlib scipy

### 6. Build the Fortran modules

If using Intel ifx: run the provided `f2py_compile.sh`.
If using gfortran (simpler): use this build script instead, saved as
`build_gfortran.sh`:

    #!/bin/bash
    set -e
    F2PY="python3 -m numpy.f2py -c"
    LIBS="-llapack -lblas"
    for mod in "static_disordered_hamiltonian datatype.f90 static_disordered_hamiltonian.f90 random_number_generator_fixedseed.f90" \
               "calculate_ll datatype.f90 calculate_ll.f90" \
               "hopping_rate datatype.f90 hopping_rate.f90" \
               "steady_state_mobility datatype.f90 steady_state_mobility.f90"; do
        set -- $mod
        name=$1; shift
        echo "Building $name..."
        $F2PY $LIBS "$@" -m $name --verbose
    done
    echo "All modules built successfully."

Make it executable and run:

    chmod +x build_gfortran.sh
    ./build_gfortran.sh

### 7. Run the code

    python polymer_intra_chain_mobility.py

Your Windows files are accessible under `/mnt/c/...` inside WSL, so you can
keep the project on your Windows drive if you prefer.

---

## OPTION B: Native Windows (No WSL)

This is more involved because Windows does not ship with a Fortran compiler
or LAPACK/BLAS. You must install them yourself.

### 1. Install Python

Download Python 3.11 or 3.12 (64-bit) from python.org. During install,
check **"Add python.exe to PATH"**.

### 2. Install a Fortran compiler and LAPACK/BLAS

**Either** install Intel oneAPI for Windows (large, but includes ifx and
MKL for LAPACK/BLAS):

- Download the Intel oneAPI Base Toolkit + HPC Toolkit for Windows.
- After install, the environment is initialized with `setvars.bat`
  (see the batch script below).

**Or** install MSYS2 + MinGW gfortran (lighter):

- Install MSYS2 from msys2.org.
- In the MSYS2 terminal:

      pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-lapack mingw-w64-x86_64-openblas

- Add `C:\msys64\mingw64\bin` to your Windows PATH.

### 3. Set up Python packages

Open Command Prompt (CMD):

    py -3.11 -m venv myenv
    myenv\Scripts\activate
    pip install numpy==1.26.4 matplotlib scipy

### 4. Build the Fortran modules

Save the batch script (below) as `build.bat` in the project folder, then in
CMD (with your venv activated) run:

    build.bat

### 5. Run the code

    python polymer_intra_chain_mobility.py

---

## Running the f2py Build Script on CMD

### Saving the file

1. Open a plain-text editor (Notepad, or better, VS Code / Notepad++).
2. Paste the batch-script contents.
3. Save with the name **`build.bat`** — in Notepad, set "Save as type" to
   **"All Files"** (not "Text Documents"), otherwise it becomes
   `build.bat.txt` and will not run.
4. Save it in the same folder as the `.f90` source files.

### Running the file

1. Open Command Prompt (search "cmd" in the Start menu).
2. Navigate to the project folder, e.g.:

       cd C:\Users\YourName\liverpool_code_project

3. Activate your Python environment:

       myenv\Scripts\activate

4. Run the build script by typing its name:

       build.bat

   (You do not type `python build.bat` — a `.bat` file runs directly.)

If it succeeds, you will see "All modules built successfully." and new
`.pyd` files (the compiled modules) will appear in the folder. These `.pyd`
files are the Windows equivalent of the Linux `.so` modules.

---

## The Windows Batch Build Script (build.bat)

Use the **intelvem** version if you installed Intel oneAPI; use the
**gfortran** version if you installed MSYS2/MinGW.

### Intel oneAPI version

    @echo off
    setlocal enabledelayedexpansion

    REM Adjust this path to your oneAPI install location
    call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

    set NPY_DISTUTILS_APPEND_FLAGS=1
    set F2PY=python -m numpy.f2py -c --fcompiler=intelvem --f77exec=ifx --f90exec=ifx
    set F90FLAGS=-check bounds -traceback -O0 -g
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

### gfortran (MSYS2/MinGW) version

    @echo off
    setlocal enabledelayedexpansion

    set F2PY=python -m numpy.f2py -c --fcompiler=gnu95
    set LIBS=-llapack -lblas

    echo Building static_disordered_hamiltonian...
    %F2PY% %LIBS% datatype.f90 static_disordered_hamiltonian.f90 random_number_generator_fixedseed.f90 -m static_disordered_hamiltonian --verbose
    if errorlevel 1 goto error

    echo Building calculate_ll...
    %F2PY% %LIBS% datatype.f90 calculate_ll.f90 -m calculate_ll --verbose
    if errorlevel 1 goto error

    echo Building hopping_rate...
    %F2PY% %LIBS% datatype.f90 hopping_rate.f90 -m hopping_rate --verbose
    if errorlevel 1 goto error

    echo Building steady_state_mobility...
    %F2PY% %LIBS% datatype.f90 steady_state_mobility.f90 -m steady_state_mobility --verbose
    if errorlevel 1 goto error

    echo All modules built successfully.
    goto end

    :error
    echo Build failed.
    exit /b 1

    :end

---

## Troubleshooting

- **"'python' is not recognized"** — Python is not on PATH. Reinstall Python
  with "Add to PATH" checked, or use `py` instead of `python`.
- **"cannot find -llapack"** — LAPACK/BLAS not installed or not on the
  linker path. On MinGW, confirm the openblas/lapack packages installed and
  that `C:\msys64\mingw64\bin` is on PATH. With Intel, use `-lmkl_rt`.
- **"No module named numpy"** — the venv is not activated, or numpy is not
  installed in it. Activate with `myenv\Scripts\activate` and
  `pip install numpy==1.26.4`.
- **ifx not found** — `setvars.bat` was not called, or the path in the
  script is wrong. Check your oneAPI install location.
- **.bat opens in a text editor instead of running** — it was saved as
  `build.bat.txt`. Rename it to `build.bat`.
- **numpy 2.x errors** — pin numpy with `pip install numpy==1.26.4`.

When in doubt, use WSL (Option A). The native Windows Fortran + LAPACK
toolchain is the main source of difficulty, and WSL avoids it entirely.
