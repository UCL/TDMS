<!-- \cond
 -->
<!-- 👆 this comment and the endcond below, tells doxygen to ignore the badges
and title at the top of README.md when building the project page (the title
would be duplicated) everything else in README.md is also the project homepage. -->

# TDMS · [![latest release](https://badgen.net/github/release/UCL/TDMS)](https://github.com/UCL/TDMS/releases)  [![license](https://badgen.net/github/license/UCL/TDMS)](https://github.com/UCL/TDMS/blob/main/LICENSE) [![Build and test](https://github.com/UCL/TDMS/actions/workflows/ci.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/ci.yml) [![MATLAB tests](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml)

> **Warning**
> This repository is a _work in progress_. The API will change without notice

<!-- \endcond -->

# Time-Domain Maxwell Solver

TDMS, the Time Domain Maxwell Solver, is a hybrid C++ and MATLAB tool for solving
Maxwell's equations to simulate light propagation through a medium. See the
[pdf documentation](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf) for
further details.

![The normed z-component of the H field incident on a cylinder](doc/assets/HzNormBanner.png)


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
where lines need to be commented in and the paths modified if cmake cannot
(1) find MATLAB, (2) find FFTW or (3) install to the default install prefix.

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

Once the executable has been compiled and installed, `tdms` should be in the `PATH`.
Check that installation worked with

```bash
$ tdms -h
```

You can invoke it directly or call it from a MATLAB script.
We recommend that beginners with MATLAB installed start with the demonstration MATLAB script.

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
-fd, --finite-difference:	Use the finite-difference solver, instead of the pseudo-spectral method.
-q:	Quiet operation. Silence all logging
-m:	Minimise output file size by not saving vertex and facet information
```

The basic workflow is with two arguments; an input file containing source fields, material composition, material properties and other simulation parameters as detailed in the documentation, and an output file name to be created.
The [`iterate_fdtd_matrix.m`](./tdms/matlab/iteratefdtd_matrix.m) script can be used to produce a valid input file from a MATLAB script and source-field functions. 

You can choose two possible solver methods: either pseudo-spectral time-domain (PSTD, the default) or finite-difference (FDTD, with option `--finite-difference`).

### Parallelism

TDMS is parallelised with [OpenMP](https://en.wikipedia.org/wiki/OpenMP). The maximum
number of threads can be set with the `OMP_NUM_THREADS` environment variable.
For example, to use 4 threads, in a bash shell, use:

```bash
$ export OMP_NUM_THREADS=4
```

Before calling the `tdms` executable.

## Citation

If you used TDMS in your research and found it helpful, please cite this work.
<!-- [zenodo/FIXME](https://zenodo.org/) -->

<!-- If you use TDMS in your work and have examples that you would like to share with other users, please get in touch with us at -->
<!-- [contact_address)[mailto:FIXME] -->
<details>
<summary>BibTEX</summary>

```bibtex
@software{tdms,
    author  = {Munro, Peter and others},
    license = {GPL-3.0},
    title   = {{TDMS - Time Domain Maxwell Solver}},
    URL     = {https://github.com/UCL/TDMS}
}
```

</details>
<details>
<summary>LaTeX</summary>

```tex
\bibitem{tdms}
P. Munro, et al \emph{TDMS - The Time-Domain Maxwell Solver}, \url{https://github.com/UCL/TDMS}.
```

</details>

## Acknowledgements

The TDMS source code was released under a GPL-3.0 License as part of a joint project between University College London's [Medical Physics and Biomedical Engineering](https://ucl.ac.uk/medphys) and [Centre for Advanced Research Computing](https://ucl.ac.uk/arc) with generous funding from the [Royal Society](https://royalsociety.org).

![medphys](doc/assets/biomedlogo.png)&nbsp;![arc](doc/assets/arclogo.png)

Development of this software has previously benefited from funding from the [Commonwealth Scholarships Commission](https://cscuk.fcdo.gov.uk/about-us/scholarships-and-fellowships/), the [Engineering and Physical Sciences Research Council](https://www.ukri.org/councils/epsrc/), and the [Australian Research Council](https://www.arc.gov.au/).

## Want to contribute?

We're grateful for bug reports, feature requests and pull requests. Please see our [contribution guidelines](https://github-pages.ucl.ac.uk/TDMS/md__c_o_n_t_r_i_b_u_t_i_n_g.html).
