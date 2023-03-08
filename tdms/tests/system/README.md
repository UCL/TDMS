# TDMS System Tests

**The system tests are currently being overhauled to provide a better workflow. This readme is subject to abrupt changes and updates whilst we perform this process.**

The system tests for `tdms` are specified through configuration files in the `config_files` directory, (creatively) named `config_XX.yaml`. The `XX` in the configuration file name should match the ID of the system tests, which are named `arc_XX`, and should also match the `test_id` field in the configuration file itself (see [below](#test-config-files)). The _reference outputs_ or _reference data_ are a collection of `.mat` files that are produced by the `tdms` executable when provided a set collection of _reference inputs_, which are also `.mat` files. These reference inputs are in turn generated via a collection of `MATLAB` scripts and functions, which build inputs such as the incident electric field and computational grid as specified in the documentation.

Each system test is named by the aforementioned `arc_XX` convention. A given system test `arc_XX` itself may consist of several _executions_ (or _runs_) of the `tdms` executable. Typically every system test has at least two runs; one for when there is no scattering object present, and one for when there is an obstacle. Other causes of additional runs might be due to the use of band-limited interpolation over cubic interpolation, for example. Each of these runs in turn has a reference input and reference output.

## Test Config Files

We use `.yaml` files to specify system tests and the parameters that must be passed into them. The first key that should be present in an input file is `test_id`, which should have a value of a string that provides the system test with a unique identifier - in our case, `config_XX.yaml`, the configuration file for the system test `arc_XX`, begins with
```yaml
test_id: 'XX'
```
Note that the quotes are necessary to ensure that the `test_id` is not interpreted as an `int`, `bool`, etc.

After the `test_id` key, all remaining keys should provide the filenames of the `.mat` input files that need to be produced for this system test to perform all of its runs. There are two ways to declare that a particular `.mat` input file needs to be regenerated; either by providing this information explicitly, or specifying another `.mat` input file to adjust.

### Providing (re)generation information

To declare an input `.mat` file `foo.mat` that needs to be generated, the configuration file should contain the key "`foo`". Within the `foo` value, we then specify additional parameters we need to know before we regenerate the input data. These parameters are detailed below:
```yaml
foo:
  # Declares that the file foo.mat needs to be (re) generated for this test to execute its runs.
  # The remainder of this block contains properties that only apply when generating, or running tests using foo.mat as an input to tdms

  # The file that is passed to iteratefdtd_matrix in run_bscan. These paths are relative to the data/input_generation/input_files directory
  input_file: input_file_name.m
  # One of fs (freespace), sph (sphere), cyl (cylinder), or sc (point-source at the origin). The shape of the scattering obstacle.
  obstacle: fs
  # The radius of obstacle in microns. For sph, the radius of the sphere. For cyl, the radius of the circular faces. For sc and fs, this option is ignored. Defaults to 15e-6 if not set.
  obstacle_radius: 5e-6
  # The refractive index of the scattering object. Not used in fs simulations. Defaults to 1.42 if not set.
  refind: 1.42
  # Whether iteratefdtd_matrix must be called in illsetup mode to setup the illumination file, prior to its call in filesetup mode. Defaults to False if not present or unpopulated
  illsetup: True
  # Whether the time-domain field needs to be computed prior to defining the scattering matrix and other material properties. Defaults to False if not present.
  calc_tdfield: True
  # One of pstd (pseudo-spectral time domain), fdtd (finite difference time domain). The time-propagation solver method to use. Defaults to fdtd if not present
  solver_method: pstd
  # One of cubic (cubic) or bli (band-limited). The interpolation method to use when extracting field values at spatial locations off the computational grid. Defaults to cubic if not present.
  interpolation: cubic

  # This key starts another block, which contains the details of all the runs within this system test that use foo.mat as their input
  runs:
    # Each key in this block defines the name of one run of tdms, using foo.mat as the input data
    run_name_1:
      # The name of the reference .mat output to compare the output of the local tdms run to.
      reference: reference.mat
      # A bool indicating whether or not tdms should use the cubic interpolation switch -c or not in this run. Defaults to False if not present.
      cubic_interpolation: True
      # A bool indicating whether or not tdms should use the fdtd solver method or not (in which case pstd is used). Defaults to False if not present.
      fdtd_solver: False
    run_name_2:
      # This is another run of tdms using the same foo.mat input
      reference: reference_2.mat
      fdtd_solver: False
```

### "Adjust"-ing another input file

The interpolation method and time-propagation solver methods used by `tdms` are controlled by two particular variables within the input `.mat` file. If one wishes to run `tdms` using the same numerical input data (incident fields, spatial grid, etc) but to toggle one of these options, it is inefficient to run the entire `run_bscan` function again. Instead, one can specify the "`adjust`" key rather than the "`input_file`" key after declaring a `.mat` input file:
```yaml
bar:
  # Declares that the file bar.mat needs to be (re) generated for this test to execute its runs.
  # The remainder of this block contains properties that only apply when generating, or running tests using bar.mat as an input to tdms

  # The adjust keyword tells the program to copy the .mat input whose name matches the value of this key. In this instance, bar.mat will copy foo.mat, plus any alterations that we specify in the remainder of the block
  adjust: foo
  # We can then overwrite the interpolation and solver_method settings that were specified in foo.mat's declaration
  # Note that if foo.mat _does not_ use the default value for (EG) `interpolation`, you _must_ explicitly specify `interpolation` as the default value here. In an adjust block, the absence of keys provides a default behaviour of "do not change" rather than "reset".
  interpolation: cubic
  solver_method: fdtd
  # We then declare runs using bar.mat as an input in the usual manner
  runs:
    bar_run_1:
      reference: bar_reference.mat
```

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

The _proposed_ workflow of a particular system test `arc_XX` is:
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

Due to licensing issues regarding running `MATLAB` on GitHub runners, we cannot utilise `matlabengine` to regenerate the reference input data during CI. (Although we are currently thinking of removing the `MATLAB` dependency which will then enable us to resolve this issue). The work-in-progress `test_regen.py` workflow can still be run locally through `pytest`, however in addition to `requirements.txt` you will also need to [install `matlabengine`](https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). See the [MathWorks page](https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html) link for detailed requirements, however you will need (at least) `Python >= 3.8` and a licensed version of `MATLAB`.

## Generating Input Data for the System Tests

The system tests rely on `.mat` input files that are presently generated through a series of MATLAB function calls and scripts. This directory contains the functionality to automatically (re)generate this input data locally, which serves two purposes:
- It means that the `.mat` _input_ files do not need to be uploaded and fetched from Zenodo each time the system tests are to be run. They can be generated locally instead. Note that the reference output files corresponding to these inputs still need to be downloaded from Zenodo.
- We are able to reliably track changes to the way we decide to handle the inputs to `tdms`, and the system tests now also provide us with a check against unexpected behaviours due to input changes.

### (Re)generation of the Data

At a glance, (re)generating the input data for a particular test case, `arc_XX`, is a three step process:
1. Determine variables, filenames, and the particular setup of `arc_XX`. For example, is an illumination file required? What are the spatial obstacles? What is the solver method? This information is stored in a `config_tc.yaml` file.
1. Call the `run_bscan.m` function (and appropriate supporting functions in the `./matlab` folder) using the information provided from the first step to produce the `.mat` input files. Each test case is also pointed to an input file (`input_file_XX.m`) which defines test-specific variables (domain size, number of period cells, material properties, etc) which are too complex to specify in a `.yaml` file.
1. Clean up the auxillary `.mat` files that are generated by this process. In particular, any `gridfiles.mat`, illumination files, or other `.mat` files that are temporarily created when generating the input `.mat` file.

### Contents of the `data/input_generation` Directory (and subdirectories)

The `run_bscan` function is inside the `bscan/` directory.

The `matlab/` directory contains functions that `run_bscan` will need to call on during the creation of the input data. This in particular includes the `iteratefdtd_matrix` function, which does the majority of the work in setting up gridfiles, illumination files, and the `.mat` inputs themselves.

The `bscan_arguments.py` file contains Python classes that handle calling the `run_bscan` function in `MATLAB` through `Python`. In time, these classes will absorb the `run_bscan` method when the `MATLAB` header dependency is removed. The `generation_data.py` file defines the epynonomous class, one instance of which has the functionality to regenerate the input data for one system test. Since the system test framework uses `pytest`, but the data generation requires `MATLAB` (for now), we use `Python` to read in and process the information that each test requires, and then call `run_bscan` with the appropriate commands from within Python.

The `regenerate_all.py` file will work through all of the `config_XX.yaml` files in the directory and regenerate the input `.mat` data corresponding to each.

The remaining `config_XX.yaml` and `input_file_XX.m` files are as mentioned in [the previous section](#regeneration-of-the-data). These contain the information about each test that Python and `run_bscan` will need to regenerate the input files.
