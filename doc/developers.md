# Developer documentation

Welcome to the developer documentation for TDMS, the Time Domain Maxwell Solver.
These pages contain the C++ API documentation as well as some other useful information for developers.

- 🐣 If you are a new user, you probably want to start with the project's [README.md](https://github.com/UCL/TDMS/blob/main/README.md) on github. There are installation and "getting started" instructions.
- 🐛 If you've spotted a bug or want to request a feature, we have a short [CONTRIBUTING.md](https://github.com/UCL/TDMS/blob/main/CONTRIBUTING.md). Then please [send an issue through github](https://github.com/UCL/TDMS/issues).
- If you want to contribute to the development (🚀!) and have read [CONTRIBUTING.md](https://github.com/UCL/TDMS/blob/main/CONTRIBUTING.md), you're in the right place.

## What?

The physics background and numerical methods used by TDMS are described in [this pdf document](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf) (the source is in [doc/latex](https://github.com/UCL/TDMS/blob/main/doc/latex)).
In terms of the code, `tdms` is a C++ binary that can be installed in your path.
We build with [CMake](https://cmake.org/) to be moderately platform-independent. We run continuous integration builds and tests on Windows, Ubuntu, and MacOS.

### Dependencies

We depend on ...

* [OpenMP](https://en.wikipedia.org/wiki/OpenMP) for parallelization.
* [fftw](https://www.fftw.org/) to calculate numerical derivatives (numerical_derivative.h)
* [spdlog](https://github.com/gabime/spdlog) for logging (this is installed automatically by a CMake [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) if not found).
* [MATLAB](https://www.mathworks.com/products/matlab.html) since we read in and write out MATLAB format files. We are actually thinking of removing this as a strict dependency ([#70](https://github.com/UCL/TDMS/issues/70)).

### Code style and other admin

We try to stick to [Google style C++](https://google.github.io/styleguide/cppguide.html) wherever practicable. <!-- And we run `clang-format` and `cppclean` linters as part of our CI. -->
There is one exception: we prefer `#pragma once` at the top of each header (because it's one line rather than three and is almost the standard anyway).

Please keep doxygen-style comments in the **headers**.
```{.cpp}
/**
 * Add two numbers together
 *
 * @param a: real number
 * @param b: real number
 * @return sum of a and b
 */
double add(double a, double b);
```
And please document the header itself with at least:
```{.cpp}
/**
 * @file filename.h
 *
 * Describe the functions/classes in here.
 */
 #pragma once
```
You can quickly check the documentation for your changes look correct with:
```{.sh}
doxygen doc/Doxygen
firefox html/index.html # or your web browser of choice
```
You should be able to find and read what you've changed.
Don't worry about doxygen for the source files (although obviously please do write helpful comments there).
We _have_ been putting doxygen comments in the [unit-test](#unit-testing) source files wherever sensible.

For Python code (e.g. in the [system tests](#system-tests)) we use [black](https://black.readthedocs.io/en/stable/) to enforce the code style.
To apply automatic code styling to staged changes in git we recommend [`pre-commit`](https://pre-commit.com/).
If you don't have it already:
```{.sh}
python -m pip install pre-commit
```

Then from the root repository directory you can add the pre-commit hooks with
```{.sh}
ls .git/hooks
pre-commit install
ls .git/hooks
```

### Compiling and debugging {#compiling}

Once you've checked the code out, compile with:
```{.sh}
cd tdms
mkdir build; cd build
cmake .. \
# -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2019b/ \
# -DFFTW_ROOT=/usr/local/fftw3/ \
# -DCMAKE_INSTALL_PREFIX=$HOME/.local/
# -DDERIVATIVE_TYPE=FD
make install
```
You may need to help CMake find MATLAB/fftw etc.

- There are two methods to choose from: pseudospectral time domain (PSTD) or finite-difference time domain (FDTD) derivatives. **This choice must be made at compile time**, to comply with MATLAB. The default is PSTD. You can select FDTD with `-DDERIVATIVE_TYPE=FD`.
- By default, build testing is turned off. You can turn it on with `-DBUILD_TESTING=ON`.
- Also by default, debug printout is off. Turn on with `-DCMAKE_BUILD_TYPE=Debug` or manually at some specific place in the code with:
```{.cpp}
#include <spdlog/spdlog.h>

// ...

spdlog::set_level(spdlog::level::debug);

// ...

spdlog::debug("Send help");
```

### Compiling on UCL's Myriad cluster
<details>
  > **Warning**
  > These instructions are a bit experimental. Please use with care (and report anything that's wrong here)!

  If you want to test changes on UCL's [Myriad](https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/) (and/or don't have MATLAB on your pesonal machine) you can try these instructions.
  Firstly, you will probably want to [forward your ssh agent](https://stackoverflow.com/questions/12257968/) for your github ssh key.
  To do this, you first need to run the following your _local_ machine:
  ```{.sh}
  ssh-add -L # check your ssh agent is running
  ssh-add /path/to/your/github/key/id_rsa
  ssh -o ForwardAgent=yes your_user@myriad.rc.ucl.ac.uk
  ```

  And once you're on Myriad:

  ```{.sh}
  git clone git@github.com:UCL/TDMS.git

  module purge
  module load beta-modules
  module load gcc-libs/9.2.0 compilers/gnu/9.2.0 xorg-utils matlab/full/r2021a/9.10 fftw/3.3.6-pl2/gnu-4.9.2 cmake/3.21.1
  cd TDMS/tdms
  mkdir build; cd build
  cmake .. \
  # -DGIT_SSH=ON
  make install
  ```

  If you get the following error (or similar)
  ```
  fatal: unable to access 'https://github.com/gabime/spdlog/': error setting certificate verify locations:
  CAfile: /etc/ssl/certs/ca-certificates.crt
  CApath: none
  ```
  it's because the MATLAB module is interfering with the SSL certificates (and we clone over https by default). This issue is known and reported. As a workaround, we've added the build option `-DGIT_SSH=ON` to switch to `git clone` over ssh instead.
</details>


## Where's the main?

The C++ `main` function is in openandorder.cpp <!-- words with a dot in them are assumed to be files so this will hyperlink to openandorder.cpp iff *that* file is also documented. --> however this only really does file I/O and setup.
The main FDTD algorithm code is in iterator.cpp <!-- won't be linked as an undocumented file doesn't exist for doxygen... this is fine, we can link to the real file in github.--> and classes included therein.

## Testing

We have two [levels of tests](https://en.wikipedia.org/wiki/Software_testing#Testing_levels): unit tests, and full system tests.

### Unit {#unit-testing}

The unit tests use [catch2](https://github.com/catchorg/Catch2/blob/devel/docs/Readme.md#top) macros. See [tests/unit](https://github.com/UCL/TDMS/blob/main/tdms/tests/unit) for good examples in the actual test code.

To write a new test, as a rough sketch you need:

```{.cpp}
/**
 * @file test_file.cpp
 * @brief Short description of the tests.
 */
#include <catch2/catch_test_macros.hpp>
#include "things_to_be_tested.h"

/**
 * @brief Detailed description of the testing.
 *
 * Maybe go into details about the test setup.
 */
TEST_CASE("Write a meaningful test case name") {
    // set up function calls or whatever
    REQUIRE_THROW(<something>)
    CHECK(<something>)
}
```
The doxygen-style comments will be included in this developer documentation.

To run the unit tests, [compile](#compiling) with `-DBUILD_TESTING=ON`. Then run `ctest` from the build directory or execute the test executable `./tdms_tests`.

It's good practice, and reassuring for your pull-request reviewers, if new C++ functionality is at covered by unit tests.

### System {#system-tests}

The full system tests are written in Python 3, and call the `tdms` executable for known inputs and compare to expected outputs.
We use [pytest](https://docs.pytest.org) and our example data is provided as zip files on [zenodo](https://zenodo.org/).

There are a few [python packages you will need](https://github.com/UCL/TDMS/blob/main/tdms/tests/requirements.txt) before running the tests so run:
```{.sh}
python -m pip install -r tdms/tests/requirements.txt
```
if you don't already have them.

When you run the tests for the first time, the example data will be downloaded to `tdms/tests/system/data` (which is [ignored by git](https://github.com/UCL/TDMS/blob/main/.gitignore)).
Subsequent runs of the test will not re-download unless you manually delete the zip file.

A good example of running the `tdms` executable for a given input and expected output is [test_arc01.py](https://github.com/UCL/TDMS/blob/main/tdms/tests/system/test_arc01.py)

You need to [compile](#compiling) `tdms`, then the system tests can be run, e.g. from the build directory:

```{.sh}
pytest ../tests/system/ -m "not fdtd_build_only"
```

A slight complication - which we hope to simplify - it is possible to build in two different configurations for two different methods: pseudospectral time domain (PSTD) or finite-difference time domain (FDTD) derivatives.
We therefore [mark](https://docs.pytest.org/en/stable/example/markers.html) some of the system tests `fdtd_build_only`.

To test the FDTD build: [recompile](#compiling) with `-DDERIVATIVE_TYPE=FD` and run:

```{.sh}
pytest ../tests/system/ -m fdtd_build_only
```
