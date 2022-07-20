import os

from pathlib import Path
from utils import work_in_zipped_dir, HDF5File, run_tdms, download_data


ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_01.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/6838866/files/arc_01.zip",
                  to=ZIP_PATH)


@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """
    Ensure that the free space PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_01/pstd_fs_input.mat", "pstd_fs_output.mat")

    reference = HDF5File("arc_01/pstd_fs_reference_output.mat")
    assert HDF5File("pstd_fs_output.mat").matches(reference)


@work_in_zipped_dir(ZIP_PATH)
def test_cyl():
    """
    Ensure that the cylinder PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_01/pstd_cyl_input.mat", "pstd_cyl_output.mat")

    reference = HDF5File("arc_01/pstd_cyl_reference_output.mat")
    assert HDF5File("pstd_cyl_output.mat").matches(reference)
