import os
from pathlib import Path

import pytest
from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_09.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/7097063/files/arc_09.zip", to=ZIP_PATH)


@pytest.mark.skip(reason="We know about this: likely a data problem.")
@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """
    Ensure that the free space PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_09/in/pstd_fs.mat", "pstd_fs_output.mat")

    reference = HDF5File("arc_09/out/pstd_fs.mat")
    assert HDF5File("pstd_fs_output.mat").matches(
        reference
    ), "The free space PSTD output with exdetintegral doesn't match the reference."


@work_in_zipped_dir(ZIP_PATH)
def test_sc():
    """
    Ensure that the (FIXME ask Peter - SC) PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_09/in/pstd_sc.mat", "pstd_sc_output.mat")

    reference = HDF5File("arc_09/out/pstd_sc.mat")
    assert HDF5File("pstd_sc_output.mat").matches(
        reference
    ), "The sc PSTD output with exdetintegral doesn't match the reference."
