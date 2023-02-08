# Input File Generation for System Tests

Tests are organised into subdirectories of the same name as the given test.
Each directory contains:
- `run_{pstd,fdtd}_bscan.m` : The "_generating_" script that generates the input `.mat` files that are to be passed to the system test(s)
- `{pstd,fdtd}_input_file.m` : A file defining the inputs to the simulation that the test case is running, to be read by `run_{pstd,fdtd}_bscan.m`
- `config.yaml` : The file that contains the run information for this test case. See `../read_config.py` for details.
- Any other auxillary files that are required for the inputs to be successfully generated.

The input files (and supplimentary `gridfiles` and `illumination` files) are not tracked by `.git` due to being binary files, however these will be placed into the particular test directory should the generating script be executed locally.

## (Re)Generating the Input Files

To generate the input files for a particular test, run the generating script in each subdirectory.
The scripts assume that they are being run in the particular test directory, so to run all at once will require changing directory between executing MATLAB scripts.
The config files will remain consistent _unless_ the generated file names are changed in the generating script, in which case these will need to be _manually_ updated.
