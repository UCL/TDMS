# Input File Generation for System Tests

Tests are organised into subdirectories of the same name as the given test.
Each directory contains:
- `run_{pstd,fdtd}_bscan.m` : The "_generating_" script that generates the input `.mat` files that are to be passed to the system test(s)
- `{pstd,fdtd}_input_file.m` : A file defining the inputs to the simulation that the test case is running, to be read by `run_{pstd,fdtd}_bscan.m`
- `config.yaml` : The file that contains the run information for this test case. See `../read_config.py` for details.
- Any other auxillary files that are required for the inputs to be successfully generated. Usually, functions in the `tdms/matlab` folder are used to generate this input data, however in certain cases it is necessary to overwrite these functions for a given test. In such cases, the "overwritten" MATLAB functions are available within the directory for the test that requires them. MATLAB's method of searching for function-files ensures that it prefers these "overwritten" functions to the "common" functions in the `tdms/matlab` directory.

## (Re)Generating the Input Files

To regenerate the `.mat` files that are to be passed into `tdms` in the system tests, run the `regenerate_input_files.sh` script from the `data_generation/` directory. This script _does not_ alter the corresponding `config.yaml` files, nor the reference data `.mat` files if they are present locally.

The input files (and supplimentary `gridfiles` and `illumination` files) are not tracked by `.git` due to being binary files, however these will be placed into each respective test directory. See the section on [regenerating the zip files](#regenerating-the-zip-files) for more information.

## Regenerating the `.zip` Files

The system test data and reference outputs are presently hosted on [`Zenodo`](https://zenodo.org/record/7440616/), available as `.zip` files that the `test_system.py` script will fetch (if they are not already present locally) when the tests run. Each `.zip` file (for each test) should contains
- The reference outputs, as `.mat` files.
- The input `.mat` files, which are those generated [from the step above](#regenerating-the-input-files).
- `config.yaml`, as outlined at the top of the page.

Having [generated the input files](#regenerating-the-input-files) and placed the appropriate reference data `.mat` files into the test directories, one can run the  `zip_test_data.sh` to create the `.zip` files that `test_system.py` is expecting, and place them directly into the expected location on the file-system if it exists. This also allows them to be easily uploaded to `Zenodo` as a new version, if changes to the codebase have necessitated a change to the reference input/outputs.
