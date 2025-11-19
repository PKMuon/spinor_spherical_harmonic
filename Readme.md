# SpinorSphericalHarmonic

## Introduction

This repository provides a C++ shared library for computing spinor spherical harmonics and associated wave functions in momentum space. 

Source:
- `spinor_spherical_harmonic.cpp`: implements the core functionality for calculating spinor spherical harmonics.
- `spinor_spherical_harmonic_example.f90`: a Fortran example program demonstrating how to complile with cpp source together.
- `spinor_spherical_harmonic_valid.f90` : a Fortran validation program to verify the correctness, be used with provided python source.
- `spinor_wave_function.cpp`: computes the wave functions in momentum space using the spinor spherical harmonics.
- `spinor_wave_function_example.f90`: a Fortran example program demonstrating how to compile with cpp source together.

Other files:
- `wavefunction_plot.ipynb`: a Jupyter notebook for visualizing and validating the computed wave functions.
- `validation.py`: a Python script to validate the results against known solutions ($$\sum_{m}Y_{+\kappa m}(\hat{p})Y^*_{+\kappa m} (\hat{p})= \frac{2j+1}{8\pi}$$).
- `convert.cpp`: converts `.dat` files to binary format so that they can be efficiently read by main programs.

## Build

requirement:

- Boost C++ Libraries

To build all source files, run:
```bash
make
```

### Optional Individual Build

To build the shared library, run:

```bash
g++ -g -O3 -shared -fPIC spinor_spherical_harmonic.cpp -o libspinor_spherical_harmonic.so
g++ -g -O3 -shared -fPIC spinor_wave_function.cpp -o libspinor_wave_function.so -lboost_math_c99
```

To build the Fortran example programs, run:

```bash
gfortran -g -O3 spinor_spherical_harmonic_example.f90 spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_example -lstdc++
gfortran -g -O3 spinor_spherical_harmonic_valid.f90 spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_valid -lstdc++
gfortran -g -O3 spinor_wave_function_example.f90 spinor_wave_function.cpp -o spinor_wave_function_example -lstdc++ -lboost_math_c99
```

To build excutable for wave function calculation, run:

```bash
g++ -g -O3 -DNEED_TEST_MAIN spinor_spherical_harmonic.cpp -o spinor_spherical_harmonic_test
g++ -g -O3 -DNEED_TEST_MAIN spinor_wave_function.cpp -o spinor_wave_function_test -lboost_math_c99
```

### Convert

To build the toolkit to convert `.dat` files to binary format, run:

```bash
g++ -std=c++17 -o convert convert.cpp
./convert input.dat output.bin
```