import os
from pathlib import Path

from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_02.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/6838977/files/arc_02.zip", to=ZIP_PATH)


@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """Ensure that the free space PSTD output matches the reference"""

    run_tdms("arc_02/pstd_fs_input.mat", "pstd_fs_output.mat")

    reference = HDF5File("arc_02/pstd_fs_reference_output.mat")
    assert HDF5File("pstd_fs_output.mat").matches(reference)


@work_in_zipped_dir(ZIP_PATH)
def test_cyl():
    """Ensure that the cylinder PSTD output matches the reference."""

    run_tdms("arc_02/pstd_cyl_input.mat", "pstd_cyl_output.mat")

    reference = HDF5File("arc_02/pstd_cyl_reference_output.mat")
    assert HDF5File("pstd_cyl_output.mat").matches(reference)
