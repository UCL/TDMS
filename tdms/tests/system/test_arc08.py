import os
from pathlib import Path

from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir

ZIP_PATH = Path(os.path.dirname(os.path.abspath(__file__)), "data", "arc_08.zip")

if not ZIP_PATH.exists():
    download_data("https://zenodo.org/record/7086087/files/arc_08.zip", to=ZIP_PATH)


@work_in_zipped_dir(ZIP_PATH)
def test_fs():
    """Ensure that the free space stady-state PSTD output matches the reference"""

    run_tdms("arc_08/in/pstd_fs_steady_fast.mat", "pstd_fs_steady_output_fast.mat")

    reference = HDF5File("arc_08/out/pstd_fs_steady_fast.mat")
    assert HDF5File("pstd_fs_steady_output_fast.mat").matches(
        reference
    ), "The free space, steady-state PSTD output doesn't match the reference."


@work_in_zipped_dir(ZIP_PATH)
def test_sph():
    """Ensure that the spherical steady-state PSTD output matches the reference."""

    run_tdms("arc_08/in/pstd_sph_steady_fast.mat", "pstd_sph_steady_output_fast.mat")

    reference = HDF5File("arc_08/out/pstd_sph_steady_fast.mat")
    assert HDF5File("pstd_sph_steady_output_fast.mat").matches(
        reference
    ), "The spherical steady-state PSTD output doesn't match the reference."
