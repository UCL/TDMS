import os
from typing import Literal

import numpy as np
import pytest

from ..bscan.composition_matrix import build_composition_matrix
from ..matlab_engine.matlab_engine import create_engine_for_testing

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))

OBSTACLES_TO_FILES = {
    "fs": os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_01.m"
    ),
    "cyl": os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_01.m"
    ),
    "sph": os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_03.m"
    ),
    "sc": os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_09.m"
    ),
}


@pytest.mark.skip(
    "MATLAB to Python indexing issue. ind2sub is odd, but the I matrix each test generates is the same, so not sure what's playing up here..."
)
@pytest.mark.parametrize(
    "obstacle, input_file_for_matlab", list(OBSTACLES_TO_FILES.items())
)
def test_build_composition_matrix(
    obstacle: Literal["fs", "cyl", "sph", "sc"], input_file_for_matlab: str
) -> None:
    radius = 15e-6

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    # Fetch the spatial grid from the matlab input file
    x, y, z, _ = engine.fdtd_bounds(input_file_for_matlab, nargout=4)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    nx = x.size
    ny = y.size
    nz = z.size

    # Compute the result according to matlab
    mat_composition_matrix = engine.composition_matrix_builder(
        input_file_for_matlab, obstacle, radius, nargout=1
    )
    mat_composition_matrix = np.array(mat_composition_matrix)
    engine.quit()

    # Compute the composition matrix according to Python
    py_composition_matrix = build_composition_matrix(x, y, z, obstacle, radius)

    if obstacle != "fs":
        # Adjust MATLAB result (in cases where the result is non-empty) for comparibility
        # - Reverse order of 1st 3 columns because of MATLAB's cloumn-major format vs numpy's row-major. Composition matrix comes out of MATLAB in the form [k, j, i, 1] when interpretted in Python, and we want [i, j, k, 1].
        # - Adjust for 1-indexing in the i, j, k columns
        # mat_composition_matrix[:,0:3] = np.flip(mat_composition_matrix[:,0:3], axis=1)
        mat_composition_matrix[:, 0:3] -= 1

        # Assert correct shape
        assert py_composition_matrix.shape == mat_composition_matrix.shape
        # Assert identical indices found
        assert np.alltrue(py_composition_matrix == mat_composition_matrix)
        # Assert hanging column of 1s is present
        assert all(py_composition_matrix[:, 3] == 1)
    else:
        # Assert we have built the empty array
        assert py_composition_matrix.size == 0

    return
