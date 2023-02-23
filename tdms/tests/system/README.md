# TDMS System Tests

**The system tests are currently being overhauled to provide a better workflow. This readme is subject to abrupt changes and updates whilst we perform this process.**

The system tests for `tdms` are specified through configuration files in the `data/input_generation/` directory, (creatively) named `config_XX.yaml`. The `XX` in the configuration file name should match the ID of the system tests, which are named `arc_XX`, and should also match the `test_id` field in the configuration file itself (see [below]()). The _reference outputs_ or _reference data_ are a collection of `.mat` files that are produced by the `tdms` executable when provided a set collection of _reference inputs_, which are also `.mat` files. These reference inputs are in turn generated via a collection of `MATLAB` scripts and functions, which build inputs such as the incident electric field and computational grid as specified in the documentation.

Each system test is named by the aforementioned `arc_XX` convention. A given system test `arc_XX` itself may consist of several _executions_ (or _runs_) of the `tdms` executable. Typically every system test has at least two runs; one for when there is no scattering object present, and one for when there is an obstacle. Other causes of additional runs might be due to the use of band-limited interpolation over cubic interpolation, for example. Each of these runs in turn has a reference input and reference output.

## Running the System Tests

The system tests can be run locally in a virtual environment that has the `requirements.txt` installed, and on a system in which `MATLAB` is installed and is on the `$PATH` under the alias `matlab`. Simply invoke `pytest` on the `tdms/tests/system` directory:
```bash
pytest  tdms/tests/system/
```

The `test_system.py` file contains the workflow that executes each system test in sequence.

**The `test_regen.py` and `tdms_testing_class.py` will be replacing the `test_system.py` and `read_config.py` files when the input overhaul is complete.**

### Workflow of a System Test

The current workflow of the system tests is:
- Fetch the reference outputs and inputs from [Zenodo](https://zenodo.org/record/7440616/files), if not already present
- Loop through each system test `arc_XX` in sequence:
    - Execute each run in `arc_XX`
    - Compare the output of each run to the corresponding reference data
    - A given run fails if the output produced by executing `tdms` differs significantly from the reference data
- `arc_XX` fails if any one of its runs fail. Otherwise it passes.

**The `test_regen.py` and `tdms_testing_class.py` will be replacing the `test_system.py` and `read_config.py` files when the input overhaul is complete.**

The proposed workflow of a particular system test `arc_XX` is:
- Locally generate the reference inputs using the `data/input_generation` functionality.
    - `arc_XX` fails if its reference input cannot be successfully generated. This indicates a failure in the scripts and/or functions in the `data/input_generation/{bscan,matlab}` directories.
- Fetch the reference outputs from [Zenodo](https://zenodo.org/record/7440616/files).
- For each run, named `run_id` in `arc_XX`:
    - Execute the call to `tdms` corresponding to `run_id`.
    - Compare the output of each run to the corresponding reference data.
    - `run_id` fails if the output produced differs significantly from the reference data.
    - Outputs produced by `run_id` are cleaned up.
- `arc_XX` fails if any one of its runs fail. Failed runs are reported by name.
- Reference inputs are cleaned up.
- `arc_XX` passes if this step is reached successfully.

## Generating Input Data for the System Tests

The system tests rely on `.mat` input files that are presently generated through a series of MATLAB function calls and scripts. This directory contains the functionality to automatically (re)generate this input data locally, which serves two purposes:
- It means that the `.mat` _input_ files do not need to be uploaded and fetched from Zenodo each time the system tests are to be run. They can be generated locally instead. Note that the reference output files corresponding to these inputs still need to be downloaded from Zenodo.
- We are able to reliably track changes to the way we decide to handle the inputs to `tdms`, and the system tests now also provide us with a check against unexpected behaviours due to input changes.

### (Re)generation of the Data

At a glance, (re)generating the input data for a particular test case, `arc_tc`, is a three step process:
1. Determine variables, filenames, and the particular setup of `arc_tc`. For example, is an illumination file required? What are the spatial obstacles? What is the solver method? This information is stored in a `config_tc.yaml` file.
1. Call the `run_bscan.m` function (and appropriate supporting functions in the `./matlab` folder) using the information provided from the first step to produce the `.mat` input files. Each test case is also pointed to an input file (`input_file_tc.m`) which defines test-specific variables (domain size, number of period cells, material properties, etc) which are too complex to specify in a `.yaml` file.
1. Clean up the auxillary `.mat` files that are generated by this process. In particular, any `gridfiles.mat`, illumination files, or other `.mat` files that are temporarily created when generating the input `.mat` file.

### Contents of the `data/input_generation` Directory (and subdirectories)

The `run_bscan` function is inside the `bscan/` directory.

The `matlab/` directory contains functions that `run_bscan` will need to call on during the creation of the input data. This in particular includes the `iteratefdtd_matrix` function, which does the majority of the work in setting up gridfiles, illumination files, and the `.mat` inputs themselves.

The `generate_test_input.py` file contains `.py` files that the system tests can invoke to regenerate the input data. Since the system test framework uses `pytest`, but the data generation requires `MATLAB` (for now), we use `Python` to read in and process the information that each test requires, and then call `run_bscan` with the appropriate commands from within Python.

The `regenerate_all.py` file will work through all of the `config_tc.yaml` files in the directory and regenerate the input `.mat` data corresponding to each.

The remaining `config_XX.yaml` and `input_file_XX.m` files are as mentioned in [the previous section](#regeneration-of-the-data). These contain the information about each test that Python and `run_bscan` will need to regenerate the input files.
