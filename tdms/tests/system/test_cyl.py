import os

from pathlib import Path
from utils import work_in_zipped_dir, HDF5File, run_tdms


THIS_DIR_PATH = os.path.dirname(os.path.abspath(__file__))


@work_in_zipped_dir(Path(THIS_DIR_PATH, "pstd_cyl.zip"))
def test_cyl():
    """
    Ensure that the cylinder PSTD output matches the reference. Checks
    that the output .mat file (with a HDF5 format) contains tensors with
    relative mean square values within numerical precision of the reference.
    """

    run_tdms("pstd_cyl/pstd_cyl_input.mat", "pstd_cyl_output.mat")

    reference = HDF5File("pstd_cyl/pstd_cyl_reference_output.mat")
    assert HDF5File("pstd_cyl_output.mat").matches(reference)
