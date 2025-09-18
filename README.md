<!--! \cond -->

# TDMS · [![latest release](https://badgen.net/github/release/UCL/TDMS)](https://github.com/UCL/TDMS/releases) [![license](https://badgen.net/badge/license/GPL-3.0/blue)](https://github.com/UCL/TDMS/blob/main/LICENSE) [![Build and test](https://github.com/UCL/TDMS/actions/workflows/ci.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/ci.yml) [![MATLAB tests](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml/badge.svg)](https://github.com/UCL/TDMS/actions/workflows/matlab_tests.yml) [![codecov](https://codecov.io/gh/UCL/TDMS/branch/main/graph/badge.svg?token=3kqP14kslL)](https://codecov.io/gh/UCL/TDMS) [![DOI](https://zenodo.org/badge/448864310.svg)](https://zenodo.org/badge/latestdoi/448864310)

<!--! @endcond -->

# Time-Domain Maxwell Solver

TDMS, the Time Domain Maxwell Solver, is a hybrid C++ and MATLAB tool for simulating light propagation through a medium by solving Maxwell's equations.
For further details about the method, please refer to the [PDF documentation](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf).

![The normed z-component of the H field incident on a cylinder](doc/assets/HzNormBanner.png)

## Getting started

We don't ship binaries at the moment, so to use TDMS, it has to be compiled.
It needs to be built against [FFTW](https://www.fftw.org/) and [MATLAB](https://www.mathworks.com/products/matlab.html), which must be downloaded and installed first.

<details>
<summary><img src="https://github.com/EgoistDeveloper/operating-system-logos/blob/master/src/24x24/UBT.png"/> Linux prerequisite setup</summary>

Assuming you don't already have them, you'll need a C++ compiler, CMake, OpenMP and FFTW.

For Debian-based distributions this should be as simple as

```{sh}
$ sudo apt install git gcc cmake libfftw3-dev libgomp1
```

<!-- TODO: add libhdf5-dev here when updating to v1.1+ -->

</details>

<details>
<summary><img src="https://github.com/EgoistDeveloper/operating-system-logos/blob/master/src/24x24/MAC.png"/> MacOS prerequisite setup</summary>

On MacOS everything can be installed using [Homebrew](https://brew.sh):

```{sh}
$ brew install cmake fftw llvm hdf5
```

</details>

<details>
<summary><img src="https://github.com/EgoistDeveloper/operating-system-logos/blob/master/src/24x24/WIN.png"/> Windows prerequisite setup</summary>

TDMS was developed - and extensively tested - on linux.
Support for Windows is quite new and experimental (please [report](https://github.com/UCL/TDMS/issues/new/choose) any issues you encounter!).

It might be more straightforward to use the [Windows subsystem for Linux (WSL2)](https://learn.microsoft.com/en-gb/windows/wsl/install), or set up an linux virtual machine.

However, TDMS _can_ be compiled natively on Windows.
This has been tested Windows 10 and 11, with PowerShell.

Assuming you don't already have them, you'll need to download and install:

<!-- when updating to version 1.1+, uncomment

* [HDF5](https://portal.hdfgroup.org/display/HDF5/HDF5+CPP+Reference+Manuals) -->

- [MATLAB](https://www.mathworks.com/products/matlab.html),
- [Visual Studio](https://visualstudio.microsoft.com/vs/community/) and **be sure to select the C++ developer kit**,
- [CMake](https://cmake.org/download/),
- and [FFTW](https://www.fftw.org/install/windows.html).

You can check that the Visual Studio compiler was installed with:

```{pwsh}
PS> MSBuild.exe -h
```

Potentially the simplest way to get FFTW is via [conda](https://anaconda.org/conda-forge/fftw):

```{pwsh}
PS> conda install -c conda-forge fftw --yes
```

<!-- TODO: add HDF5 👆 here when updating to v1.1+ probably also via conda is the easiest -->

You'll need to ensure the paths to FFTW and MATLAB (the locations of `fftw3.dll` and `libmex.dll` respectively) are in the `env:Path`.

These can be found, e.g. by

```{pwsh}
PS> conda list fftw # assuming you installed via conda
PS> which.exe MATLAB
```

Which should return something like `C:\Program Files (x86)\MATLAB\R20XXx\bin\matlab` and maybe `C:\ProgramData\envs\base\bin`.
If you downloaded FFTW and created `fftw3.dll` with `lib.exe`, you just need to know where you saved it.

You can append the paths:

```{pwsh}
PS> $env:Path += ";C:\Program Files (x86)\MATLAB\R20XXx\bin\;C:\ < wherever fftw3.dll is >"
```

Which will help Windows locate `.dll` files later.
For all following instructions, you'll have to substitute our mentions of `tdms` with `tdms.exe` and `$` is used to denote a command prompt which, in PowerShell, would look like `PS>`

<details>
<summary>Even more Windows troubleshooting</summary>

We've seen that in a fresh PowerShell window, Windows does not remember the location of the `.dll` files, so you may have to re-add them to the path, or copy them into the directory where TDMS was installed.

TDMS typically installs to `"C:\Program Files (x86)\tdms\bin\tdms.exe"`.

</details>

</details>

---

You'll need to download and install [MATLAB](https://www.mathworks.com/products/matlab.html), and take note where the headers are installed.

</details>

## Installation

To compile and install, follow these steps:

```{sh}
$ git clone git@github.com:UCL/TDMS.git
$ cd TDMS
$ git checkout v1.0.2 # the stable version
$ mkdir build; cd build
$ cmake ../tdms \
$ # -DMatlab_ROOT_DIR=/usr/local/MATLAB/R20XXx/ \
$ # -DFFTW_ROOT=/usr/local/fftw3/ \
$ # -DCMAKE_INSTALL_PREFIX=$HOME/.local/
$ cmake --build . --target install --config Release
```

If CMake cannot find MATLAB, FFTW, or install to the default installation prefix, uncomment the relevant line(s) and modify the path(s) accordingly.

<!-- TODO: add HDF5 when updating to v1.1+ -->

You can check that `tdms` was installed correctly and is in your `PATH` by running:

```{sh}
$ tdms --help
$ tdms --version
```

in a new terminal.

## How to run

You can run TDMS either directly or from a MATLAB script.
For beginners, we recommend starting with the demonstration MATLAB script, which you can find in the `examples/arc_01` directory.
Move into this directory, launch MATLAB, and run the MATLAB script [`run_pstd_bscan.m`](https://github.com/UCL/TDMS/blob/main/examples/arc_01/run_pstd_bscan.m).
This script will generate the input to TDMS, run TDMS, and display sample output.
There are comments explaining what it is doing, so you can follow along with what is being setup and created at each stage.
We have also commented the input file [`arc_01_example_input.m`](https://github.com/UCL/TDMS/blob/main/examples/arc_01/arc_01_example_input.m) that this script passes to `iteratefdtd_matrix.m`.

<details>
<summary>Troubleshooting</summary>

We've seen that launching MATLAB on MacOS via the launcher (cmd + space) may not preserve the system `PATH`.

```
command not found: tdms
```

Assuming `tdms --help` works in a new terminal, try launching MATLAB _from_ that terminal.

```{sh}
$ tdms --help
$ /Applications/MATLAB_<version>.app/bin/matlab
```

The MATLAB example scripts should then find `tdms`.
If you still have problems, you can try hard-coding the full path to `tdms` into the MATLAB script.

In a terminal run

```{sh}
$ which tdms
```

Copy the full path (something like `/usr/local/bin/tdms`) into [`run_pstd_bscan.m`](https://github.com/UCL/TDMS/blob/main/examples/arc_01/run_pstd_bscan.m), replacing the `'tdms'` text in the calls to the `system()` function.

</details>

### MATLAB file version

In order to be readable by TDMS, files need to be saved in .mat (MATLAB file) version 7.3 or newer.
This can be done by passing '-v7.3' to MATLAB's save command as the final argument, for example:

```
save('my_input_file.mat', 'my_input_var_1', 'my_input_var_2', ..., 'my_input_var_N', '-v7.3');
```

<details>
<summary>Or set v7.3 as the default in MATLAB's settings.</summary>

![](doc/assets/matlab-file-settings.png)

</details>

### On the command line

If you want to run TDMS standalone at the command line, the basic operation is with two arguments: an input file containing simulation parameters, and an output file name.
You can choose between two solver methods: **finite-difference** or **pseudo-spectral**, as well as two interpolation methods: **cubic** or **bandlimited**.
These options can be selected by setting the corresponding flag variables in the input file.
When `tdms` reads the input, it will verify if the input file contains a dataset that matches the names of these flags.

There are two flags available for configuration in the input file.

<details>
<summary> `use_pstd` </summary>
- If not provided, or provided as `false`, then the default timestepping method of finite-differences (FDTD) will be used.
- If present and set to `true`, then `tdms` will use the pseudo-spectral (PSTD) method when performing simulation timesteps.
</details>
<details>
<summary> `use_bli` </summary>
- If not provided, or provided as `false`, then the default interpolation method of cubic interpolation will be used to obtain field values at the centres of Yee cells.
- If present and set to `true`, then `tdms` will use bandlimited interpolation (BLI) when obtaining field values at Yee cell centres.

\note Typically bandlimited interpolation is superior to cubic interpolation when the extent of the Yee cell is of approximately the same order as, but slightly less than, one-sixth of the shortest wavelength of interest.
Otherwise, cubic interpolation typically enjoys superior accuracy.

</details>

TDMS is parallelised with [OpenMP](https://en.wikipedia.org/wiki/OpenMP).
You can set the maximum number of threads using the `OMP_NUM_THREADS` environment variable before calling the TDMS executable.

```{sh}
$ export OMP_NUM_THREADS=4 # for example
```

## Citation

If you used TDMS in your research and found it helpful, please cite this work: [10.5281/zenodo.7950604](https://doi.org/10.5281/zenodo.7950604).

<!-- If you use TDMS in your work and have examples that you would like to share with other users, please get in touch with us at -->
<!-- [contact_address)[mailto:FIXME] -->
<details>
<summary>BibTEX</summary>

```bibtex
@software{tdms,
    author       = {Munro, Peter and others},
    license      = {GPL-3.0},
    title        = {{TDMS - The Time-Domain Maxwell Solver}},
    URL          = {https://github.com/UCL/TDMS},
    publisher    = {Zenodo},
    doi          = {10.5281/zenodo.7950603}
}
```

</details>
<details>
<summary>LaTeX</summary>

```tex
\bibitem{tdms}
P. Munro, et al \emph{TDMS - The Time-Domain Maxwell Solver}, \url{https://github.com/UCL/TDMS}, \href{https://doi.org/10.5281/zenodo.7950603}{10.5281/zenodo.7950603}.
```

</details>

## Acknowledgements

The TDMS source code was released under a GPL-3.0 License as part of a joint project between University College London's [Medical Physics and Biomedical Engineering](https://ucl.ac.uk/medphys) and [Centre for Advanced Research Computing](https://ucl.ac.uk/arc) with generous funding from the [Royal Society](https://royalsociety.org).

![medphys](doc/assets/biomedlogo.png)&nbsp;![arc](doc/assets/arclogo.png)

Development of this software has previously benefited from funding from the [Commonwealth Scholarships Commission](https://cscuk.fcdo.gov.uk/about-us/scholarships-and-fellowships/), the [Engineering and Physical Sciences Research Council](https://www.ukri.org/councils/epsrc/), and the [Australian Research Council](https://www.arc.gov.au/).

## Want to contribute?

We're grateful for bug reports, feature requests, and pull requests. Please see our [contribution guidelines](https://github-pages.ucl.ac.uk/TDMS/md__c_o_n_t_r_i_b_u_t_i_n_g.html) (we also have some [developer documentation](https://github-pages.ucl.ac.uk/TDMS/md_doc_developers.html)).
