import os
from glob import glob
from pathlib import Path

import pytest
from data.input_generation.generate_test_input import (
    generate_test_input as regenerate_test,
)
from run_system_test import run_system_test
from utils import download_data

# Location of this file, which is where the tests are running from
LOCATION_OF_THIS_FILE = Path(os.path.abspath(os.path.dirname(__file__)))
# This will determine whether or not we want to retain the regenerated input .mat files (if say, we are planning a new Zenodo upload). Recommended FALSE on CLI, TRUE locally if you're doing the update
PRESERVE_FLAG = True

# Location of the config files that contain test information
path_to_config_files = Path(LOCATION_OF_THIS_FILE, "config_files")
# The data/generation directory
path_to_input_generation = Path(LOCATION_OF_THIS_FILE, "data", "input_generation")

# Dataset is stored at https://zenodo.org/record/7440616/
# ccaegra@ucl.ac.uk (William Graham, @willGraham01) has access.
ZENODO_URL = "https://zenodo.org/record/7440616/files"
# directory in which to store the downloaded zip files
ZIP_DESTINATION = Path(os.path.dirname(os.path.abspath(__file__)), "data")

# Find all config_XX.yaml files and hence infer all test cases that need to be run
TEST_IDS = sorted(glob(str(path_to_config_files / "config_*.yaml")))
for i, id in enumerate(TEST_IDS):
    # remove everything bar the ID from what we found
    id = id.removeprefix(str(path_to_config_files / "config_"))
    id = id.removesuffix(".yaml")
    TEST_IDS[i] = id
# The number of tests will be useful to know later
N_TESTS = len(TEST_IDS)


def workflow(test_id: str, preserve_inputs: bool = PRESERVE_FLAG) -> None:
    """Performs all tdms runs contained in the system test arc_{test_id}. Assumes that reference data is present in the appropriate .zip folder within the /data directory.

    The workflow steps are as follows:
    1. Regenerate the input data for arc_{test_id}
    2. Perform each run of TDMS as specified by arc_{test_id}, comparing the output to the previously obtained references
    3. Perform clean-up on generated outputs and inputs, so they are not saved to the cache and accidentally re-used
    """
    # Fetch config file location
    config_file_path = path_to_config_files / f"config_{test_id}.yaml"

    # Regenerate the input data for arc_{test_id}
    # Error on not successful (hence fail test)
    print(f"Regenerating input data...", end="", flush=True)
    regenerate_test(config_file_path)
    print(f"Done")

    # Perform each run of TDMS as specified by the config file
    run_success = run_system_test(config_file_path)

    # TEAR-DOWN: remove the regenerated inputs and outputs from the data/input_generation/arc_{test_id} folder.
    # To be safe, we can just remove all .mat files from this directory the the subdirectories, since these should be the only system-test TDMS artefacts
    if not preserve_inputs:
        input_data_dump = str(path_to_input_generation / f"arc_{test_id}/*.mat")
        mat_artefacts = glob(input_data_dump, recursive=True)
        for mat_file in mat_artefacts:
            os.remove(mat_file)

    # Although we should have check-ed whether each run was a pass/fail, we can also assert that all runs need to pass here to report failures
    failed_run_names = []
    for run_id, passed in run_success.items():
        if not passed:
            failed_run_names.append(run_id)
    assert all(
        run_success
    ), f"arc_{test_id} : Some runs were unsuccessful ({len(run_success)-sum(run_success)}/{len(run_success)}) :\n {failed_run_names}"
    return


# Run each system test that is currently saved in the repository, as identified through the config files
# "MATLAB doesn't run on GH runners. If you want to check that regenerating the input data still allows the tests to pass, run this locally and remove the skip mark."
@pytest.mark.parametrize("test_id", TEST_IDS)
def test_system(test_id) -> None:
    """Runs the system test arc_{test_id}, including fetching missing reference data from Zenodo.

    Wraps the workflow() method, which actually does the bulk of the testing.
    """
    string_to_indicate_test_start = f"Running arc_{test_id}"
    test_header_padding = "=" * len(string_to_indicate_test_start)
    print(f"\n{test_header_padding}")
    print(string_to_indicate_test_start)
    print(test_header_padding)

    # the reference OUTPUT data should be at this location
    ZIP_PATH = ZIP_DESTINATION / f"arc_{test_id}.zip"
    # download data if not present
    if not ZIP_PATH.exists():
        url = ZENODO_URL + f"/arc_{test_id}.zip"
        print(f"Downloading from {url}")
        download_data(url, to=ZIP_PATH)
    else:
        print(f"Using cache in {ZIP_PATH}")

    # Run the test workflow, for this test
    workflow(test_id)
    # End of system test
    return
