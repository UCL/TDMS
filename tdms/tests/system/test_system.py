from pathlib import Path

import pytest
from utils import compare_output, download_data, work_in_zipped_dir

ZENODO_URL = "https://zenodo.org/record/"
TEST_URLS = {
    "01": ZENODO_URL + "6838866/files/arc_01.zip",
    "02": ZENODO_URL + "6838977/files/arc_02.zip",
    "03": ZENODO_URL + "6839280/files/arc_03.zip",
}


@pytest.mark.parametrize("number", TEST_URLS.keys())
@pytest.mark.parametrize("run_type", ["fs", "cyl", "sph"])
def test_system(number: str, run_type: str):
    """
    Run the system tests. For each of the test data URLs defined above this:
    - downloads the test data
    - runs TDMS with the input file
    - compares the output to the reference output file
    """
    ZIP_PATH = (Path(__file__).parent / "data" / f"arc_{number}.zip").resolve()

    if not ZIP_PATH.exists():
        download_data(TEST_URLS[number], to=ZIP_PATH)

    # Need to define a new compare function to make sure we do the comparison
    # in the zipped directory
    @work_in_zipped_dir(ZIP_PATH)
    def compare(number, run_type):
        input_fname = f"arc_{number}/pstd_{run_type}_input.mat"
        output_fname = f"arc_{number}/pstd_{run_type}_reference_output.mat"

        if not Path(input_fname).exists():
            pytest.skip(
                f"No input file for run type '{run_type}' in test number {number}"
            )

        compare_output(input_fname, output_fname)

    compare(number, run_type)
