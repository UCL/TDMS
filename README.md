[![Linux tests](https://github.com/UCL/TDMS/actions/workflows/linux_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/linux_tests.yml)
[![Windows tests](https://github.com/UCL/TDMS/actions/workflows/windows_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/windows_tests.yml)
[![MATLAB tests](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml)
[![MacOS tests](https://github.com/UCL/TDMS/actions/workflows/macos_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/macos_tests.yml)
[![doc](https://img.shields.io/badge/PDF-latest-orange.svg?style=flat)](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf)
# TDMS

> **Warning**
> This repository is a _work in progress_. The API will change without notice

***
## Introduction

TDMS (Time Domain Maxwell Solver) is a hybrid C++ and MATLAB for solving
Maxwell's equations to simulate light propagation through a medium. See the
[pdf documentation](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf) for
further details.


***
## Compilation

TDMS requires building against [FFTW](https://www.fftw.org/) and
[MATLAB](https://www.mathworks.com/products/matlab.html), thus both need to be
downloaded and installed prior to compiling TDMS. Install with

```bash
cd tdms
mkdir build; cd build
cmake .. \
# -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2019b/ \
# -DFFTW_ROOT=/usr/local/fftw3/ \
# -DCMAKE_INSTALL_PREFIX=$HOME/.local/ \
# -DBUILD_TESTING=ON
make install
```
where lines need to be commented in and the paths modified if cmake cannot
(1) find MATLAB, (2) find FFTW or (3) install to the default install prefix.

- <details>
    <summary>Mac specific instructions</summary>

    To compile on a Mac an x86 compiler with libraries for OpenMP are required,
    which can be installed using [brew](https://brew.sh/) with `brew install llvm`
    then (optionally) set the following cmake arguments

    ```
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


***

## Usage

Once the executable has been compiled and installed, `tdms` should be in the `PATH`.
Check that installation worked with

```bash
$ tdms -h
```

You can invoke it directly or call it from a MATLAB script.
We recommend that beginners with MATLAB installed start with the demonstration MATLAB script.

#### To run the demonstration code

Move into directory [`examples/arc_01`](./examples/arc_01/),
launch MATLAB and run the MATLAB script:

[`run_pstd_bscan.m`](./examples/arc_01/run_pstd_bscan.m)

This script will generate the input to the executable, run the executable and
display sample output.

#### To run standalone

You can also run `tdms` from the command line...

```bash
$ tdms --help
Usage:
tdms [options] infile outfile
tdms [options] infile gridfile outfile
Options:
-h:	Display this help message
-fd, --finite-difference:	Use the finite-difference solver, instead of the pseudo-spectral method.
-q:	Quiet operation. Silence all logging
-m:	Minimise output file size by not saving vertex and facet information
```

The basic workflow is with two arguments, an input file as specified by [`iterate_fdtd_matrix.m`](./tdms/matlab/iteratefdtd_matrix.m), and an output file name to be created.

You can choose two possible solver methods: either pseudo-spectral time-domain (PSTD, the default) or finite-difference (FDTD, with option `--finite-difference`).

#### Parallelism

TDMS is parallelised with [OpenMP](https://en.wikipedia.org/wiki/OpenMP). The maximum
number of threads can be set with the `OMP_NUM_THREADS` environment variable.
For example, to use 4 threads, in a bash shell, use:

```bash
$ export OMP_NUM_THREADS=4
```

Before calling the `tdms` executable.

## Contributing

You are welcome to report new issues or submit pull requests.  For more information about how to contribute, please see the [`CONTRIBUTING.md`](./CONTRIBUTING.md) file and we have developer documentation on our [`gh-pages` site](https://github-pages.ucl.ac.uk/TDMS).
