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
