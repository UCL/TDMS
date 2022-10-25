import os
from pathlib import Path

from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_13.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/7233118/files/arc_13.zip", to=ZIP_PATH)


@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """
    Ensure that the free space FDTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("--finite-difference", "arc_13/in/fdtd_fs.mat", "fdtd_fs_output.mat")

    reference = HDF5File("arc_13/out/fdtd_fs.mat")
    assert HDF5File("fdtd_fs_output.mat").matches(
        reference
    ), "The free space FDTD output doesn't match the reference."


@work_in_zipped_dir(ZIP_PATH)
def test_cyl():
    """
    Ensure that the cylinder FDTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("--finite-difference", "arc_13/in/fdtd_cyl.mat", "fdtd_cyl_output.mat")

    reference = HDF5File("arc_13/out/fdtd_cyl.mat")
    assert HDF5File("fdtd_cyl_output.mat").matches(
        reference
    ), "The spherical FDTD output doesn't match the reference."
