import os
from pathlib import Path

import pytest
from pytest_check import check
from read_config import YAMLTestConfig
from utils import HDF5File, download_data, work_in_zipped_dir

ZENODO_URL = "https://zenodo.org/record/"
# directory in which to store the downloaded zip files
ZIP_DESTINATION = Path(os.path.dirname(os.path.abspath(__file__)), "data")

# all test cases and where to aquire their data
TEST_URLS = {
    "01": ZENODO_URL + "6838866/files/arc_01.zip"  # this is still the old data!
}


def workflow() -> None:
    config_for_tests = YAMLTestConfig()
    system_tests_in_this_folder = config_for_tests.run_list

    for system_test in system_tests_in_this_folder:
        # run test
        system_test.run()
        # load reference file
        reference_file = HDF5File(system_test.get_reference_file_name())
        # load output file
        output_file = HDF5File(system_test.get_output_file_name())
        # assert file contents match
        # but continue to run the other tests in this zip file, even if we fail
        with check:
            test_passed, comparison_information = output_file.matches(
                reference_file, return_message=True
            )
            assert (
                test_passed
            ), f"In {config_for_tests.test_id} -> {system_test.run_id}: {comparison_information}"
    return


# run this for each key, value pair in the TEST_URLS (IE each test_id and zenodo url)
@pytest.mark.parametrize("test_id, url", list(TEST_URLS.items()))
def test_system(test_id, url) -> None:
    print(f"Running {test_id}", end=" | ")

    # the data should be at this location
    ZIP_PATH = ZIP_DESTINATION / f"new_arc_{test_id}.zip"
    # download data if not present
    if not ZIP_PATH.exists():
        print(f"Downloading from {url}")
        download_data(url, to=ZIP_PATH)
    else:
        print(f"Using cache in {ZIP_PATH}")

    # run the workflow in this zip folder
    work_in_zipped_dir(ZIP_PATH)(workflow)()
    return
