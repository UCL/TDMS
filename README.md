[![Run tests](https://github.com/UCL/TDMS/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/UCL/TDMS/actions/workflows/tests.yml)
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
# -DCMAKE_INSTALL_PREFIX=$HOME/.local/
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

#### To run the demonstration code

Once the executable has been compiled, move into directory _examples/arc_01_,
launch Matlab and run the Matlab script:

_run_pstd_bscan.m_

This script will generate the input to the executable, run the executable and 
display sample output.


#### Parallelism

TDMS is parallelised with [OpenMP](https://en.wikipedia.org/wiki/OpenMP). The maximum 
number of threads can be set with the `OMP_NUM_THREADS` environment variable. 
For example, to use 4 threads, in a bash shell, use:

```bash
export OMP_NUM_THREADS=4
```


## Contributing

You are welcome to report new issues or submit pull requests.  For more information about how to contribute, please see the [`CONTRIBUTING.md`](./CONTRIBUTING.md) file.
