# Developer documentation

Welcome to the developer documentation for TDMS, the Time Domain Maxwell Solver.
These pages contain the C++ API documentation as well as some other useful information for developers.

## Developers' introduction to TDMS

The physics background and numerical methods used by TDMS are described in [this pdf document](https://github.com/UCL/TDMS/blob/gh-doc/masterdoc.pdf) (the source is in [doc/latex](https://github.com/UCL/TDMS/blob/main/doc/latex)).
In terms of the code, `tdms` is a C++ binary that can be installed in the `PATH`.
We build with [CMake](https://cmake.org/) to be moderately platform-independent.
We run continuous integration builds and tests on Windows, Ubuntu, and MacOS.

## Dependencies

We depend on ...

* [OpenMP](https://en.wikipedia.org/wiki/OpenMP) for parallelization.
* [fftw](https://www.fftw.org/) to calculate numerical derivatives (numerical_derivative.h)
* [spdlog](https://github.com/gabime/spdlog) for logging (this is installed automatically by a CMake [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) if not found).
* [MATLAB](https://www.mathworks.com/products/matlab.html) since we read in and write out MATLAB format files. We are actually thinking of removing this as a strict dependency ([#70](https://github.com/UCL/TDMS/issues/70)).

### Versions

We follow [semantic versioning](https://semver.org): major, minor, patch.
Major versions may contain changes to the API and UI, minor versions contain new features but preserve backwards-compatibility of UI, patches are bugfixes.

So (for example) version v1.2.0 is major version 1, minor version 2, and the zeroth patch. v1.2.0 will be backwards-compatible with v1.1.99.

**Mechanics**: Version control is via `git tag`s on the `main` branch and propagated to the C++ source [code via CMake](https://github.com/UCL/TDMS/blob/main/tdms/cmake/version_from_git.cmake).
Release notes are kept in the [GitHub release page](https://github.com/UCL/TDMS/releases).

As well as official releases, a placeholder "version" is assigned to every commit on every development branch, or can be supplied at compile-time (which takes the highest precedence).

```{.sh}
$ cmake .. <usual cmake options> -DTDMS_VERSION=this_version_string_will_be_highest_priority
$ make install
$ tdms --version
TDMS version: this_version_string_will_be_highest_priority
```

If on a development branch then the placeholder version is `<branchname>_<tag>` if tagged, or `<branchname>_<commit hash>` if not tagged. This is useful when debugging problems introduced between commits (with the assistance of `git bisect`).

```{.sh}
$ git checkout myfeaturebranch
$ cmake ..
$ make install
$ tdms --version
TDMS version: myfeaturebranch_a341ab
```

### Code style and other admin {#code-style-and-doxygen}

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


### Linters and automatic code style {#pre-commit}

For C++ we've found [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) works well.
For Python code (e.g. in the [system tests](#system-tests)) we use [black](https://black.readthedocs.io/en/stable/).

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


## Compiling and debugging {#compiling}

Once you've checked the code out, compile with:
```{.sh}
mkdir build; cd build
cmake ../tdms \
# -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2019b/ \
# -DFFTW_ROOT=/usr/local/fftw3/ \
# -DCMAKE_INSTALL_PREFIX=$HOME/.local/ \
# -DBUILD_TESTING=ON \
# -DCODE_COVERAGE=ON \
# -DCMAKE_BUILD_TYPE=Debug
make install
```
You may need to help CMake find MATLAB/fftw etc.

- By default, build testing is turned off. You can turn it on with `-DBUILD_TESTING=ON`.
- By default, code coverage is disabled. You can enable it with `-DCODE_COVERAGE=ON`.
   - [Then follow some extra steps](#coverage).
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
<summary>Myriad instructions</summary>

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

The C++ `main` function is in `main.cpp` however the main algorithm code is in the `execute()` method of the `SimulationManager` class, in `execute_simulation.cpp`.

Broadly speaking, the `main` function - and thus a call to (or execution of) `tdms` -  can be broken down into the following steps:
1. Parse the input arguments
2. [Prepare the simulation](#initialisation-step-2) by creating an instance of `SimulationManager`
3. [Run the _main loop_](#running-the-main-loop-step-3). The _main loop_ is meant to mean "the (implementation of) the time-stepping algorithm which TDMS uses", and is contained in the `execute()` method of the `SimulationManager` class.
4. [Perform post-loop processing](#post-processing-step-4) on the outputs, having completed the main loop (extraction of phasors, computing derived quantities, etc).
5. [Write the outputs](#write-out-and-tear-down-step-5) and tear-down memory.

The `SimulationManager` class governs steps 2 through 5.

## Code organisation of the TDMS algorithm

<!-- TODO We definitely need a diagram here. -->

An instance of `SimulationManager` essentially supports one `tdms` simulation, however functionality is divided across its member variables. These objects themselves group variables by functionality, purpose, and necessary scope. The list below provides a description of these objects (with parentheses denoting the member name in the `SimulationManager` instance):
- `ObjectsFromInfile` (`inputs`): Handles conversion of any input arrays to native `C++` objects, and any variables that are derived directly from the input file supplied to `tdms` on the command line.
- `OutputMatrices` (`outputs`): Handles data storage for the outputs and the process of writing the output file.
- `PSTDVariables` (`PSTD`) and `FDTDBootstrapper` (`FDTD`): Handle variables that are only required for the pseudo-spectral and finite-difference solver methods, respectively.
- `LoopTimers` (`timers`): Controls timing (and logging execution times) of the main loop and subprocesses therein.

The role of `SimulationManager` is to handle how these objects interact in the main loop.

### 2. Initialisation

`SimulationManager` instances begin by initialising the `inputs` attribute. Based on this attribute, they set up `PSTD`/`FDTD` variables and prepare the `outputs` object for writing. Some members of the `outputs` objects are used as the buffer for the field and phasor data whilst the main loop is running, with the final state of this buffer being the output values. Other attributes of the `outputs` object that require information _from_ the main loop are prepared (dimensions are determined, etc) at this stage but not set.

### 3. Running the Main Loop

Once initialised, the `execute()` method can be called. This will create an instance of `LoopVariables` that is scoped to this function; which in turn sets up the variables that are required in the main loop but not needed elsewhere. This avoids bloating the `execute()` method with a large number of variable declarations and setup, as well as simplifying tear-down at the end of the loop. `timers` track the speed of execution of the main loop, and report to the logging.

### 4. Post-Processing

Upon completing `execute()` successfully, the `post_loop_processing()` method of `SimulationManager` writes the data to the `outputs` attributes that were not used as buffers during the main loop.

### 5. Write out and Tear Down

From here, the `write_outputs()` method exports the data in `outputs` to the desired file.

At this point, the instance of `SimulationManager` can be allowed to go out of scope. In practice, `main()` terminates here and the destructor handles the tear-down of `MATLAB` memory blocks.

## Testing {#testing}

We have two [levels of tests](https://en.wikipedia.org/wiki/Software_testing#Testing_levels): unit tests, and full system tests.

### Unit {#unit-testing}

The unit tests use [catch2](https://github.com/catchorg/Catch2/blob/devel/docs/Readme.md#top) macros. See [tests/unit](https://github.com/UCL/TDMS/blob/main/tdms/tests/unit) for good examples in the actual test code.

To write a new test, as a rough sketch you need:

```{.cpp}
/**
 * @file test_file.cpp
 * @brief Short description of the tests.
 */
#include "things_to_be_tested.h"

#include <catch2/catch_test_macros.hpp>

/**
 * @brief Detailed description of the testing.
 *
 * Maybe go into details about the test setup.
 */
TEST_CASE("Write a meaningful test case name") {
    // set up function calls or whatever
    REQUIRE_THROW(<something>);
    CHECK(<something>);
}
```
The doxygen-style comments will be included in this developer documentation.

To run the unit tests, [compile](#compiling) with `-DBUILD_TESTING=ON`. Then run `ctest` from the build directory or execute the test executable `./tdms_tests`.

It's good practice, and reassuring for your pull-request reviewers, if new C++ functionality is at covered by unit tests.

### Test coverage {#coverage}

If you want to check the coverage of the unit tests, then build with `-DCODE_COVERAGE=ON`.
Then when you run `ctest` or the `tdms_tests` executable, [GCDA coverage files](https://gcc.gnu.org/onlinedocs/gcc/Gcov-Data-Files.html) are generated.
They are not really for humans to parse but can be summarised with other tools.
In our build system we use [`lcov`](https://github.com/linux-test-project/lcov), a frontend for `gcov` reports.
If you don't have `lcov`, you can install it with aptitude or homebrew (`sudo apt install lcov` or `brew install lcov`).

To get a coverage report in the terminal: execute the unit tests and then confirm you generated .gcda files...
```{.sh}
./tdms_tests
find . --name "*.gcda"
```
they're probably somewhere like `build/CMakeFiles/tdms.dir/src`.

Then generate the coverage report with
```{.sh}
lcov --capture --directory ./CMakeFiles/tdms_tests.dir/src --output-file coverage.info
lcov --list coverage.info
```
You'll notice it also counts the dependencies.
You can filter these (`lcov --remove`) if you don't want them to appear in the statistics:
```{.sh}
lcov --remove coverage.info '/usr/*' --output-file coverage.info
lcov --remove coverage.info '*/_deps/*' --output-file coverage.info
```

But you don't need to worry about this too much.
When you open a pull-request, codecov will report.
And you can browse the [codecov TDMS pages online](https://codecov.io/gh/UCL/TDMS).

### MATLAB Units {#matlab-unit-tests}

To run the unit tests, go to the `tdms/tests/matlab` directory and use the command:

```bash
matlab -batch run_tests
```

In the `utils` folder, you'll find small MATLAB functions to edit the `data/pstd_input_file{2,3}D.m` files during the testing phase.

The `iteratefdtd_matrix` function has different ways to generate source terms from an input file.
However, some combinations of variables provided can be incompatible, resulting in errors or contradictory information.
The unit test `test_iteratefdtd_matrix_source_combinations` goes through each possible combination of source-related inputs, checking for errors or valid input combinations.
The combinations and expected results are listed in the table below.

\note `TD-field` is also known as `exi/eyi` in some docstrings.

\note `usecd` was the legacy name for the variable that controlled which solver method to use in the timestepping algorithm. Previous convention was that `usecd = 1` (or not present) resulted in the use of FDTD. This has since been superceeded by the `use_pstd` flag which is `true` when PSTD is to be used, and FDTD will be used otherwise (such as when this flag is not present or set explicitly to `false`).

| TD-field | Using FDTD | compactsource | efname    | hfname    | Raises error? | Error info |
|:--------:|:----------:|:-------------:|:---------:|:---------:|:-------------:|:----------:|
| 1        | 1          | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 1          | 1             | 1         | 0         | 0             |  |
| 1        | 1          | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 1          | 1             | 0         | 0         | 0             |  |
| 1        | 1          | 0             | 1         | 1         | 0             |  |
| 1        | 1          | 0             | 1         | 0         | 1             | If not compact source, both efname and hfname must be specified |
| 1        | 1          | 0             | 0         | 1         | 1             | If not compact source, both efname and hfname must be specified |
| 1        | 1          | 0             | 0         | 0         | 0             |  |
| 1        | 0          | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 0          | 1             | 1         | 0         | 0             |  |
| 1        | 0          | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 0          | 1             | 0         | 0         | 0             |  |
| 1        | 0          | 0             | 1         | 1         | 1             | Cannot use FDTD & not compactsource |
| 1        | 0          | 0             | 1         | 0         | 1             | Cannot use FDTD & not compactsource |
| 1        | 0          | 0             | 0         | 1         | 1             | Cannot use FDTD & not compactsource |
| 1        | 0          | 0             | 0         | 0         | 1             | Cannot use FDTD & not compactsource |
| 0        | 1          | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 1          | 1             | 1         | 0         | 0             |  |
| 0        | 1          | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 1          | 1             | 0         | 0         | 1             | Must specify efname if compact source and TD-field not specified |
| 0        | 1          | 0             | 1         | 1         | 0             |  |
| 0        | 1          | 0             | 1         | 0         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 1          | 0             | 0         | 1         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 1          | 0             | 0         | 0         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 0          | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 0          | 1             | 1         | 0         | 0             |  |
| 0        | 0          | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 0          | 1             | 0         | 0         | 1             | Must specify efname if compact source and TD-field not specified |
| 0        | 0          | 0             | 1         | 1         | 1             | Cannot use FDTD & not compactsource |
| 0        | 0          | 0             | 1         | 0         | 1             | Cannot use FDTD & not compactsource |
| 0        | 0          | 0             | 0         | 1         | 1             | Cannot use FDTD & not compactsource |
| 0        | 0          | 0             | 0         | 0         | 1             | Cannot use FDTD & not compactsource |

### System {#system-tests}

The full system tests are written in Python 3, and call the `tdms` executable for known inputs and compare to expected outputs.
We use [pytest](https://docs.pytest.org) and our example data is provided as zip files on [zenodo](https://zenodo.org/).

There are a few [python packages you will need](https://github.com/UCL/TDMS/blob/main/tdms/tests/requirements.txt) before you are able to run the tests, which can be installed by executing:
```{.sh}
python -m pip install -r tdms/tests/requirements.txt
```
if you don't already have them.
You'll then need to [compile](#compiling) `tdms` with `-DBUILD_TESTING=ON`.
Once compiled, the system tests can be run by invoking `pytest` and pointing it to the `tdms/tests/system` directory.
For example, from the build directory:
```{.sh}
$ pwd
/path/to/repository/TDMS/tdms/build
$ pytest ../tests/system/
```
The [`test_system.py`](https://github.com/UCL/TDMS/blob/main/tdms/tests/system/test_system.py) script runs each system test in sequence.

**[`test_regen.py`](https://github.com/UCL/TDMS/blob/main/tdms/tests/system/test_regen.py) and [`tdms_testing_class.py`](https://github.com/UCL/TDMS/blob/main/tdms/tests/system/tdms_testing_class.py.py) will replace `test_system.py` and `read_config.py` when the [input overhaul](https://github.com/UCL/TDMS/issues/70) is complete.**

When you run the tests for the first time, test data is downloaded to `tdms/tests/system/data` (and will be [ignored by git](https://github.com/UCL/TDMS/blob/main/.gitignore)).
These reference input files contain arrays for: the incident electric field, the computational grid, etc. which are needed by the simulation, and have been generated by a trusted version of the relevant MATLAB scripts.
Subsequent runs of the tests will not re-download unless you manually delete the zip file(s).

The system tests for `tdms` are configured with yaml files in the `data/input_generation/` directory.
They are named `config_XX.yaml` where `XX` matches the ID of the system test, which themselves are named `arc_XX` by historical convention.
These config files define a number of variables, parameters, and filenames that will be required when running the system tests - for a detailed breakdown of their syntax, see the annotated `config_01.yaml` file.
This should also match the `test_id` field in the configuration file itself.
The _reference outputs_ or _reference data_ are a collection of `.mat` files, produced from the _reference inputs_ by a trusted version of the `tdms` executable.
We test for regression using these reference files.

A given system test typically has two calls to the `tdms` executable; one for when there is no scattering object present, and one for when there is some obstacle.
More than two runs in a test might be due to the use of band-limited interpolation over cubic interpolation.
Each call to the executable has a reference input and reference output.
In the scripts, a given execution is called by `tests.utils.run_tdms` which wraps a [subprocess.Popen](https://docs.python.org/3/library/subprocess.html#subprocess.Popen).

#### Workflow of a System Test

The workflow of a particular system test `arc_XX` is:
- Locally generate the reference inputs using functionality in `data/input_generation/generate_test_data.py`.
    - `arc_XX` fails if its reference input cannot be successfully generated.
- Fetch the reference outputs from [Zenodo](https://zenodo.org/record/7440616/files).
- For each run, named `run_id` in `arc_XX`:
    - Execute the call to `tdms` corresponding to `run_id`.
    - Compare the output of each run to the corresponding reference data.
    - `run_id` fails if the output produced differs significantly from the reference data.
    - Outputs produced by `run_id` are cleaned up.
- `arc_XX` fails if any one of its runs fail. Failed runs are reported by name.
- Reference inputs are cleaned up.
- `arc_XX` passes if this step is reached successfully.

Due to [licensing issues regarding running `MATLAB` on GitHub runners](https://github.com/matlab-actions/setup-matlab/issues/13), we cannot use `matlabengine` to regenerate the reference input data during CI. (Although we are currently thinking of removing the `MATLAB` dependency which will then enable us to resolve this issue). The work-in-progress `test_regen.py` workflow can still be run locally through `pytest`, however in addition to `requirements.txt` you will also need to [install `matlabengine`](https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). See the [MathWorks page](https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html) link for detailed requirements. You will need (at least) `Python >= 3.8` and a licensed version of `MATLAB`.

### Generating Input Data for the System Tests

The system tests rely on `.mat` input files that are generated through a series of MATLAB function calls and scripts. This directory contains the functionality to automatically (re)generate this input data, which serves two purposes:
-  The `.mat` _input_ files do not need to be uploaded and fetched from Zenodo each time the system tests are run. They can be generated locally instead.
   - Note that the reference output files corresponding to these inputs still need to be downloaded from Zenodo.
- We track changes to the way we handle inputs to `tdms`, and the system tests. Ensuring we test against unexpected behaviour due to input changes.

#### (Re)generation of the Data

(Re)generating the input data for a particular test case, `arc_XX`, is a three-step process:
1. Determine variables, filenames, and the particular setup of `arc_XX`. This information is stored in the corresponding `config_XX.yaml` file. For example, is an illumination file required? What are the spatial obstacles? What is the solver method?
1. Call the `run_bscan.m` function (and sub-functions in `./matlab`) using the information in `config_XX.yaml` to produce the `.mat` input files. Each test case requires an input file (`input_file_XX.m`) which defines test-specific variables (domain size, number of period cells, material properties, etc) which are too complex to specify in a `.yaml` file.
1. Clean up the auxillary `.mat` files that are generated by this process. In particular, any `gridfiles.mat`, illumination files, or other `.mat` files that are temporarily created when generating the input `.mat` file.

#### Contents of the `data/input_generation` Directory (and subdirectories)

The `run_bscan` function is inside the `bscan/` directory.

The `matlab/` directory contains functions that `run_bscan` will need to call on during the creation of the input data. This in particular includes the `iteratefdtd_matrix` function, which does the majority of the work in setting up gridfiles, illumination files, and the `.mat` inputs themselves.

The `generate_test_input.py` file contains `.py` files that the system tests can invoke to regenerate the input data. Since the system test framework uses `pytest`, but the data generation requires `MATLAB` (for now), we use `Python` to read in and process the information that each test requires, and then call `run_bscan` with the appropriate commands from within Python.

The `regenerate_all.py` file will work through all of the `config_XX.yaml` files in the directory and regenerate the input `.mat` data corresponding to each.

The remaining `config_XX.yaml` and `input_file_XX.m` files are as mentioned in [the previous section](#regeneration-of-the-data). These contain the information about each test that Python and `run_bscan` will need to regenerate the input files.
