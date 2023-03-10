<!-- \cond
 -->
<!-- 👆 this comment and the endcond below, tells doxygen to ignore the badges
and title at the top of README.md when building the project page (the title
would be duplicated) everything else in README.md is also the project homepage. -->

[![Linux tests](https://github.com/UCL/TDMS/actions/workflows/linux_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/linux_tests.yml)
[![Windows tests](https://github.com/UCL/TDMS/actions/workflows/windows_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/windows_tests.yml)
[![MATLAB tests](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml)
[![MacOS tests](https://github.com/UCL/TDMS/actions/workflows/macos_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/macos_tests.yml)
[![doc](https://img.shields.io/badge/PDF-latest-orange.svg?style=flat)](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf)
# TDMS

> **Warning**
> This repository is a _work in progress_. The API will change without notice

<!-- \endcond -->

# Time-Domain Maxwell Solver

TDMS, the Time Domain Maxwell Solver, is a hybrid C++ and MATLAB tool for
solving Maxwell's equations to simulate light propagation through a medium. See
the [pdf documentation](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf)
for further details.


## Compilation

TDMS needs to be built against [FFTW](https://www.fftw.org/) and
[MATLAB](https://www.mathworks.com/products/matlab.html), thus both need to be
downloaded and installed prior to compiling TDMS. Install with

```bash
$ cd tdms
$ mkdir build; cd build
$ cmake .. \
# -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2019b/ \
# -DFFTW_ROOT=/usr/local/fftw3/ \
# -DCMAKE_INSTALL_PREFIX=$HOME/.local/ \
# -DBUILD_TESTING=ON
$ make install
```
where lines need to be commented in and the paths modified if cmake cannot (1)
find MATLAB, (2) find FFTW or (3) install to the default install prefix.

<details>
<summary>Mac specific instructions</summary>

To compile on a Mac an x86 compiler with libraries for OpenMP are required,
which can be installed using [brew](https://brew.sh/) with `brew install llvm`
then (optionally) set the following cmake arguments

```{sh}
-DCMAKE_CXX_COMPILER=/Users/username/.local/homebrew/opt/llvm/bin/clang++
-DOMP_ROOT=/Users/username/.local/homebrew/opt/llvm/
-DCXX_ROOT=/Users/username/.local/homebrew/opt/llvm
```

On an ARM Mac install the x86 version of brew with
```bash
arch -x86_64 zsh
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
arch -x86_64 /usr/local/bin/brew install llvm
```
</details>


## Usage

Once the executable has been compiled and installed, `tdms` should be in the
`PATH`.  Check that installation worked with

```bash
$ tdms -h
```

You can invoke it directly or call it from a MATLAB script.  We recommend that
beginners with MATLAB installed start with the demonstration MATLAB script.

### To run the demonstration code

Move into directory [`examples/arc_01`](./examples/arc_01/),
launch MATLAB and run the MATLAB script:

[`run_pstd_bscan.m`](./examples/arc_01/run_pstd_bscan.m)

This script will generate the input to the executable, run the executable and
display sample output.

### To run standalone

You can also run `tdms` from the command line...

```bash
$ tdms --help
Usage:
tdms [options] infile outfile
tdms [options] infile gridfile outfile
Options:
-h:	Display this help message
-q:	Quiet operation. Silence all logging
-m:	Minimise output file size by not saving vertex and facet information
```

The basic workflow is with two arguments, an input file as specified by
[`iterate_fdtd_matrix.m`](./tdms/matlab/iteratefdtd_matrix.m), and an output
file name to be created.

### Executable / simulation options

`tdms` has a couple of options that change the algorithms that will be used at various points in the code.

#### Solver Methods

You can choose two possible solver methods to conduct the time-propagation portion of the simulation: either finite-difference time domain (FDTD) or psuedo-spectral time-domain (PSTD).
The method used is governed by the `usecd` variable in the `infile` passed to `tdms`;
- `usecd = 1` will result in FDTD being used,
- `usecd = 0` will result in PSTD being used.
If `usecd` is not provided in the `infile`, FDTD will be used.

#### Interpolation Methods

You may choose between two possible methods for interpolating field values away from the computational grid, which are used when extracting field values across a user-defined surface or at user-defined vertices.
Your choices are cubic interpolation, which uses four datapoints during interpolation and fits a cubic polynomial through them, and band-limited interpolation [FIXME] with 8 datapoints.
Band-limited interpolation provides superior accuracy over cubic interpolation when the longest wavelength of interest is of the same order, but slightly less than, the dimensions of the Yee cell over 6.
The interpolation method used is controlled by the `intmethod` variable in the `infile` passed to `tdms`;
- `intmethod = 1` will result in cubic interpolation,
- `intmethod = 2` will result in band-limited interpolation.
If `intmethod` is not provided in the `infile`, cubic interpolation will be used.

### Parallelism

TDMS is parallelised with [OpenMP](https://en.wikipedia.org/wiki/OpenMP). The
maximum number of threads can be set with the `OMP_NUM_THREADS` environment
variable.  For example, to use 4 threads, in a bash shell, use:

```bash
$ export OMP_NUM_THREADS=4
```

Before calling the `tdms` executable.

## Want to contribute?

We're very grateful for bug reports, feature requests and pull requests. Please
see our [contribution guidelines](./CONTRIBUTING.md).
