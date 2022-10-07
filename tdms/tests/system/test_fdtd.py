import os
from pathlib import Path

import pytest
from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(
    os.path.dirname(os.path.abspath(__file__)), "data", "arc_example_fdtd.zip"
)

if not ZIP_PATH.exists():
    download_data(
        "https://zenodo.org/record/7148198/files/arc_example_fdtd.zip", to=ZIP_PATH
    )


@pytest.mark.fdtd_build_only
@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """
    Ensure that the free space FDTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.

    Note: At the moment (due to MATLAB) this test only works if tdms is compiled with
    """

    # misleading filename here
    run_tdms("arc_example_fdtd/in_pstd_fs.mat", "fdtd_fs_output.mat")

    reference = HDF5File("arc_example_fdtd/out_fdtd_fs.mat")
    assert HDF5File("fdtd_fs_output.mat").matches(reference)


@pytest.mark.fdtd_build_only
@work_in_zipped_dir(ZIP_PATH)
def test_cyl():
    """
    Ensure that the cylinder FDTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    # misleading filename here
    run_tdms("arc_example_fdtd/in_pstd_cyl.mat", "fdtd_cyl_output.mat")

    reference = HDF5File("arc_example_fdtd/out_fdtd_cyl.mat")
    assert HDF5File("fdtd_cyl_output.mat").matches(reference)
