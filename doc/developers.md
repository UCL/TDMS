# Developer documentation

Welcome to the developer documentation for TDMS, the Time Domain Maxwell Solver.
These pages contain the C++ API [doxygen](https://doxygen.nl) documentation as well as some other useful information for developers.

- 🐣 If you are a new user, you probably want to start with the project's [README.md](https://github.com/UCL/TDMS/blob/main/README.md) on github. There are installation and "getting started" instructions.
- 🐛 If you've spotted a bug or want to request a feature, we have a short [CONTRIBUTING.md](https://github.com/UCL/TDMS/blob/main/CONTRIBUTING.md). Then please [send an issue through github](https://github.com/UCL/TDMS/issues).
- If you want to contribute to the development (🚀!) and have read [CONTRIBUTING.md](https://github.com/UCL/TDMS/blob/main/CONTRIBUTING.md), you're in the right place.

## What?

The physics and method used by TDMS is described in [this pdf document](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf) (the source is in [doc/latex](https://github.com/UCL/TDMS/blob/main/doc/latex)).
In terms of the code, `tdms` is a C++ binary that can be installed in your path.
We build with [CMake](https://cmake.org/) to be moderately platform-independent. We run CI builds and tests in Windows, Ubuntu, and MacOS.

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

### Compiling and debugging

- By default, build testing is turned off. You can turn it on with `-DBUILD_TESTING=ON`.
- Also by default, debug printout is off. Turn on with `-DCMAKE_BUILD_TYPE=Debug` or manually at some specific place in the code with:
```{.cpp}
spdlog::set_level(spdlog::level::debug);

// ...

spdlog::debug("Send help");
```

## Where's the main ?

The C++ `main` function is in openandorder.cpp <!-- words with a dot in them are assumed to be files so this will hyperlink to openandorder.cpp iff *that* file is also documented. --> however this only really does file I/O and setup. 
The main FDTD algorithm code is in iterator.cpp <!-- won't be linked as an undocumented file doesn't exist for doxygen... this is fine, we can link to the real file in github.--> and classes included therein.
