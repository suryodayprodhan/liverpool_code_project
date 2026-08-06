#!/bin/bash
set -e

source /opt/intel/oneapi/setvars.sh

export NPY_DISTUTILS_APPEND_FLAGS=1
export NPY_NUMPY_SITE_CFG=~/.numpy-site.cfg

F2PY="python3.11 -m numpy.f2py -c --fcompiler=intelem --f77exec=ifx --f90exec=ifx"
F90FLAGS="-O3 -xHost"
LIBS="-llapack -lblas"

echo "Building static_disordered_hamiltonian..."
$F2PY $LIBS \
    datatype.f90 \
    static_disordered_hamiltonian.f90 \
    random_number_generator_fixedseed.f90 \
    -m static_disordered_hamiltonian \
    --f90flags="$F90FLAGS" \
    --verbose

echo "Building calculate_ll..."
$F2PY $LIBS \
    datatype.f90 \
    calculate_ll.f90 \
    -m calculate_ll \
    --f90flags="$F90FLAGS" \
    --verbose

echo "Building hopping_rate..."
$F2PY $LIBS \
    datatype.f90 \
    hopping_rate.f90 \
    -m hopping_rate \
    --f90flags="$F90FLAGS" \
    --verbose

echo "Building steady_state_mobility..."
$F2PY $LIBS \
    datatype.f90 \
    steady_state_mobility.f90 \
    -m steady_state_mobility \
    --f90flags="$F90FLAGS" \
    --verbose

echo "All modules built successfully."
