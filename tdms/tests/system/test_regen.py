import os
import sys
from glob import glob
from pathlib import Path
from typing import Union

# Location of this file, which is where the tests are running from
LOCATION_OF_THIS_FILE = Path(os.path.abspath(os.path.dirname(__file__)))

# Hack the 1st - add the data/generation directory to the path so we can import from it
path_to_input_generation = Path(LOCATION_OF_THIS_FILE, "data", "input_generation")
sys.path.insert(0, str(path_to_input_generation))

import pytest
from regenerate_all import regenerate_test
from utils import download_data

# Dataset is stored at https://zenodo.org/record/7440616/
# ccaegra@ucl.ac.uk (William Graham, @willGraham01) has access.
ZENODO_URL = "https://zenodo.org/record/7440616/files"
# directory in which to store the downloaded zip files
ZIP_DESTINATION = Path(os.path.dirname(os.path.abspath(__file__)), "data")

# Find all config_XX.yaml files and hence infer all test cases that need to be run
TEST_IDS = glob(str(Path(path_to_input_generation, "config_*.yaml")))
for i, id in enumerate(TEST_IDS):
    # remove everything bar the ID from what we found
    id = id.removeprefix(str(Path(path_to_input_generation, "config_")))
    id = id.removesuffix(".yaml")
    TEST_IDS[i] = id
# The number of tests will be useful to know later
N_TESTS = len(TEST_IDS)


def workflow(test_id: str, ZIP_PATH: Union[Path, str]) -> None:
    """Performs all tdms runs contained in the system test arc_{test_id}. Assumes that reference data is present in the appropriate .zip folder within the /data directory.

    The workflow steps are as follows:
    1. Regenerate the input data for arc_{test_id}
    2. Perform each run of TDMS as specified by arc_{test_id}, comparing the output to the previously obtained references
    3. Perform clean-up on generated outputs, so they are not saved to the cache and accidentally re-used
    """
    # Fetch config file location before beginning
    config_file = Path(path_to_input_generation, f"config_{test_id}.yaml")

    # Regenerate the input data for arc_{test_id}, fail test if unsuccessful
    regeneration_success = regenerate_test(config_file)
    assert regeneration_success == 0, f"Data regeneration for arc_{test_id} failed!"

    # Perform each run of TDMS as specified by the config file - note: no longer in zipped dir so need to rework my old class :(
    # config_for_tests = YAMLTestConfig()
    # system_tests_in_this_folder = config_for_tests.run_list

    # for system_test in system_tests_in_this_folder:
    #     # run test
    #     system_test.run()
    #     # load reference file
    #     reference_file = HDF5File(system_test.get_reference_file_name())
    #     # load output file
    #     output_file = HDF5File(system_test.get_output_file_name())
    #     # assert file contents match
    #     # but continue to run the other tests in this zip file, even if we fail
    #     with check:
    #         test_passed, comparison_information = output_file.matches(
    #             reference_file, return_message=True
    #         )
    #         assert (
    #             test_passed
    #         ), f"In {config_for_tests.test_id} -> {system_test.run_id}: {comparison_information}"

    # Perform clean-up on generated outputs
    return


# Run each system test that is currently saved in the repository, as identified through the config files
@pytest.mark.parametrize("test_id", TEST_IDS)
def test_system(test_id) -> None:
    """Runs the system test arc_{test_id}, including fetching missing reference data from Zenodo.

    Wraps the workflow() method, which actually does the bulk of the testing.
    """
    print(f"\nRunning {test_id}", end=" | ")

    # the reference OUTPUT data should be at this location
    ZIP_PATH = ZIP_DESTINATION / f"arc_{test_id}.zip"
    # download data if not present
    if not ZIP_PATH.exists():
        url = ZENODO_URL + f"/arc_{test_id}.zip"
        print(f"Downloading from {url}")
        download_data(url, to=ZIP_PATH)
    else:
        print(f"Using cache in {ZIP_PATH}")

    # regenerate the input data for this particular test
    workflow(test_id, ZIP_PATH)

    # TEAR-DOWN: remove the regenerated inputs from the data/input_generation/arc_{test_id} folder.
    # To be safe, we can just remove all .mat files from this directory the the subdirectories, since these should be the only system-test TDMS artefacts
    input_data_dump = str(path_to_input_generation + f"/arc_{test_id}/*.mat")
    mat_artefacts = glob(input_data_dump, recursive=True)
    for mat_file in mat_artefacts:
        os.remove(mat_file)
    return


# I'm only here so Will can press F5 in VSCode to test this file, rather than switching windows to invoke pytest :) I will disappear soon
test_system("01")
