import glob
import os
from pathlib import Path

import pytest
from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZENODO_URL = "https://zenodo.org/record/"
"""
To add a new test, upload the test files to Zenodo and add a new
entry to this dictionary. The number on the LHS must be the same
as the number in the zip filename.

The zip file must contain:
- a single input file named `pstd_input_file.m`
- one or more input data files named `pstd_{run_type}_input.mat`.
  run_type can be any string.
- corresponding reference output data files named `pstd_{run_type}_reference_output.mat`.
  For every input file there must be a corresponding output file with the same run_type in the filename.

run_type is one of ["fs", "cyl", "sph"].
"""
TEST_URLS = {
    "01": ZENODO_URL + "6838866/files/arc_01.zip",
    "02": ZENODO_URL + "6838977/files/arc_02.zip",
    "03": ZENODO_URL + "6839280/files/arc_03.zip",
    # 08 is currently failing for known reasons
    # "08": ZENODO_URL + "7086087/files/arc_08.zip",
    "09": ZENODO_URL + "7097063/files/arc_09.zip",
    "10": ZENODO_URL + "7096040/files/arc_10.zip",
}


@pytest.mark.parametrize("number", TEST_URLS.keys())
def test_system(number: str):
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
    def compare(number):
        pattern = "pstd_*_input.mat"
        input_fnames = glob.glob(pattern)
        assert (
            len(input_fnames) > 0
        ), f"No input files found in {os.getcwd()}. Make sure input files match pattern {pattern}."
        for input_fname in input_fnames:
            output_fname = input_fname.replace("input", "reference_output")
            if not Path(output_fname).exists():
                raise RuntimeError(
                    f"Could not find expected reference file: {output_fname}"
                )
            compare_output(input_fname, output_fname)

    compare(number)


def compare_output(input_filename: os.PathLike, reference_output_filename: os.PathLike):
    """
    Run TDMS using `input_filename`, then compare the output to the data
    saved in `reference_output_filename`.

    Checks that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """
    output_filename = "pstd_output.mat"
    run_tdms(input_filename, output_filename)

    reference = HDF5File(reference_output_filename)
    output = HDF5File(output_filename)
    try:
        output.assert_matches(reference)
    finally:
        os.remove(output_filename)
