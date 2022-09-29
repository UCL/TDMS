import os
from pathlib import Path

from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_10.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/7096040/files/arc_10.zip", to=ZIP_PATH)


@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """
    Ensure that the free space PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_10/in/pstd_fs.mat", "pstd_fs_output.mat")

    reference = HDF5File("arc_10/out/pstd_fs.mat")
    assert HDF5File("pstd_fs_output.mat").matches(
        reference
    ), "The free space PSTD output with exi present doesn't match the reference."


@work_in_zipped_dir(ZIP_PATH)
def test_cyl():
    """
    Ensure that the cylinder PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("arc_10/in/pstd_cyl.mat", "pstd_cyl_output.mat")

    reference = HDF5File("arc_10/out/pstd_cyl.mat")
    assert HDF5File("pstd_cyl_output.mat").matches(
        reference
    ), "The spherical PSTD output with exi present doesn't match the reference."
