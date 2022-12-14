import os
from pathlib import Path

from read_config import YAMLTestConfig
from utils import download_data, work_in_zipped_dir

ZENODO_URL = "https://zenodo.org/record/"
# directory in which to store the downloaded zip files
ZIP_DESTINATION = Path(os.path.dirname(os.path.abspath(__file__)), "data")

# all test cases and where to aquire their data
TEST_URLS = {
    "01": ZENODO_URL + "6838866/files/arc_01.zip"  # this is still the old data!
}

# download data zip files, or use existing ones if already present
TEST_ZIP_FOLDERS = {}
for test_id, url in TEST_URLS.items():
    ZIP_PATH = ZIP_DESTINATION / f"new_arc_{test_id}.zip"
    # download data if not present
    if not ZIP_PATH.exists():
        print(f"Downloading {test_id} from {url}")
        download_data(url, to=ZIP_PATH)
    # add entry to the dictionary of tests
    TEST_ZIP_FOLDERS[test_id] = ZIP_PATH


def test_system() -> None:
    # for each zip folder containing tests...
    for test_id, ZIP in TEST_ZIP_FOLDERS.items():
        print(f"Now running {test_id} in {ZIP}")

        # this is the workflow for the test in the zip folder
        @work_in_zipped_dir(ZIP)
        def run_test_in_zip():
            config_for_tests = YAMLTestConfig()
            system_tests_in_this_folder = config_for_tests.generate_test_list()

            for system_test in system_tests_in_this_folder:
                system_test.run()
            return

        # now run said workflow
        run_test_in_zip()
    return
