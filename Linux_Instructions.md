# Running the Polymer Intra-Chain Mobility Code on Linux

This code has two parts:

1. **Fortran modules** (`static_disordered_hamiltonian`, `calculate_ll`,
   `hopping_rate`, `steady_state_mobility`) compiled with f2py into
   Python-importable modules.
2. **Python scripts** (`polymer_intra_chain_mobility.py`, `diffusivity.py`,
   `mobility_sweep.py`, etc.) that import those compiled modules.

The code is built with the **Intel oneAPI toolkit** (the `ifx` Fortran
compiler and Intel MKL for LAPACK/BLAS), which gives optimized, high-
performance numerical routines. The provided `f2py_compile.sh` is written for
this toolchain.

These instructions assume a Debian/Ubuntu-based Linux distribution. For other
distributions, replace `apt` with your package manager (`dnf`, `yum`,
`pacman`, `zypper`, etc.), and use Intel's corresponding repository.

---

## 1. System Build Tools

Install the basic compilers, Python, and the venv module:

    sudo apt update
    sudo apt install -y build-essential python3 python3-pip python3-venv python3-dev

`build-essential` provides a C compiler and `make`, which f2py needs in
addition to the Intel Fortran compiler.

---

## 2. Intel oneAPI Toolkit

The `f2py_compile.sh` script uses the Intel Fortran compiler `ifx` and Intel
MKL for LAPACK/BLAS. Install oneAPI from Intel's APT repository:

    # Add Intel's GPG key
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
        | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

    # Add the oneAPI repository
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
        | sudo tee /etc/apt/sources.list.d/oneAPI.list

    sudo apt update

    # Install the Base Toolkit (MKL) and HPC Toolkit (ifx Fortran compiler)
    sudo apt install -y intel-basekit intel-hpckit

This is a large download (several GB). Once installed, the compiler and
libraries live under `/opt/intel/oneapi/`.

### Sourcing the environment

You must source the oneAPI environment before building or running. This puts
`ifx`, MKL, and the required runtime libraries on your paths:

    source /opt/intel/oneapi/setvars.sh

To avoid typing this every session, add it to your `~/.bashrc`:

    echo 'source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1' >> ~/.bashrc

Verify the compiler is found:

    which ifx
    ifx --version

Note: recent oneAPI releases ship `ifx` (not the older `ifort`). The provided
build script already uses `--f77exec=ifx --f90exec=ifx` to handle this. If
`which ifx` returns nothing, the environment has not been sourced.

---

## 3. Python Environment (without conda)

Create an isolated virtual environment with `venv` and activate it:

    python3 -m venv gnn
    source gnn/bin/activate

Your prompt should now show `(gnn)`. To leave the environment later, run
`deactivate`.

If your system Python is older than 3.11, install a newer one first:

    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt update
    sudo apt install -y python3.12 python3.12-venv python3.12-dev
    python3.12 -m venv gnn
    source gnn/bin/activate

---

## 4. Python Packages

With the environment activated, install the required packages and versions:

    pip install --upgrade pip
    pip install numpy==1.26.4 scipy matplotlib 
    pip install meson ninja

Version notes:

- **numpy must be 1.x** (pinned to 1.26.4). The Fortran build via
  `numpy.f2py` uses `numpy.distutils`, which was removed in numpy 2.x, and
  compiled modules can crash under a numpy 2.x runtime. Do not upgrade numpy.
- **scipy** is required by `diffusivity.py` (sparse matrix support).
- **matplotlib** is required by `mobility_sweep.py` for the plots.

---

## 5. Build the Fortran Modules

Make sure the oneAPI environment is sourced and the Python environment is
activated first (so f2py uses the right numpy):

    source /opt/intel/oneapi/setvars.sh      # if not already sourced
    source gnn/bin/activate                  # if not already active

Then run the provided build script:

    chmod +x f2py_compile.sh
    ./f2py_compile.sh

The script builds all four modules and prints
"All modules built successfully." on completion. Compiled `.so` files appear
in the current directory - these are the importable modules.

### Production vs debug flags

For production runs (much faster), edit the `F90FLAGS` line in
`f2py_compile.sh` from the debug flags:

    F90FLAGS="-check bounds -traceback -O0 -g"

to optimized flags:

    F90FLAGS="-O3 -xHost"

`-O0` disables optimization and `-check bounds` checks every array access at
runtime; switching to `-O3 -xHost` often gives a large speedup and takes full
advantage of the Intel compiler's optimizations for your CPU. Keep the debug
flags only when actually debugging.

---

## 6. Run the Code

With the environment active and oneAPI sourced:

    python polymer_intra_chain_mobility.py

This reads `DEFAULT_INPUT`. To use a custom input file:

    python polymer_intra_chain_mobility.py my_input.txt

To run the parameter sweep (from inside the TestExample folder):

    cd TestExample
    python mobility_vs_disorder.py

Because the sweep launches the main program as a subprocess, it inherits the
current environment - so make sure oneAPI is sourced and the venv is active
in the shell you launch it from.

---

## 7. Quick Reference: Full Setup From Scratch

    # System tools
    sudo apt update
    sudo apt install -y build-essential python3 python3-pip python3-venv python3-dev

    # Intel oneAPI
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

---

## Troubleshooting

- **`No module named numpy`** - the venv is not activated, or numpy is not
  installed in it. Run `source gnn/bin/activate` then reinstall.
- **`ImportError: libifport.so.5: cannot open shared object file`** - the
  Intel runtime libraries are not on the library path. Source oneAPI:
  `source /opt/intel/oneapi/setvars.sh`. If it persists, add the compiler lib
  directory manually, e.g.
  `export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/latest/lib:$LD_LIBRARY_PATH`.
- **`CompilerNotFound: intelem: f90 nor f77`** - `ifx` was not found. Confirm
  oneAPI is sourced and `which ifx` returns a path. The build script uses
  `--f77exec=ifx --f90exec=ifx` to point f2py at it.
- **numpy 2.x errors / distutils missing** - numpy was upgraded past 1.x.
  Reinstall the pinned version: `pip install numpy==1.26.4`.
- **`setvars.sh: already run` warning** - harmless; it just means the
  environment was already sourced this session.
- **Permission denied running the script** - make it executable:
  `chmod +x f2py_compile.sh`.
