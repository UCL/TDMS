import os

import numpy as np
import pytest
from pytest_check import check

from ..matlab_engine.matlab_engine import create_engine_for_testing
from ..misc_functions.gauss_legendre import gauss_legendre
from ..misc_functions.yee_cell_position import yee_position

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))


def test_yee_position() -> None:
    """Test yee_position against it's MATLAB counterpart, yeeposition."""
    # MATLABEngine needs float dtype otherwise attempts integer multiplication for some reason...
    i_int = np.arange(0, 8, step=1)
    i_float = np.arange(0, 8, step=1, dtype=float)
    j_int = np.arange(2, 7, step=1)
    j_float = np.arange(2, 7, step=1, dtype=float)
    k_int = np.arange(1, 9, step=1)
    k_float = np.arange(1, 9, step=1, dtype=float)
    delta = {
        "x": 5.0e-2,
        "y": 1.0e-2,
        "z": 0.75e-2,
    }

    for component in ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]:
        with check:
            engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
            x_mat, y_mat, z_mat = engine.yeeposition(
                i_float, j_float, k_float, delta, component, nargout=3
            )
            engine.quit()

            x_mat = np.array(x_mat)
            y_mat = np.array(y_mat)
            z_mat = np.array(z_mat)
            x_py, y_py, z_py = yee_position(i_int, j_int, k_int, delta, component)

            assert np.all(
                x_py == x_mat
            ), f"Error: x co-ords for {component}-component don't match"
            assert np.all(
                y_py == y_mat
            ), f"Error: y co-ords for {component}-component don't match"
            assert np.all(
                z_py == z_mat
            ), f"Error: z co-ords for {component}-component don't match"

    return


@pytest.mark.parametrize("n_samples, a, b", [(19, -1.0, 1.0), (20, 0.0, 1.0)])
def test_gauss_legendre(n_samples: int, a: float, b: float) -> None:
    """gauss_legendre should be equivalent to it's MATLAB counterpart, gauss_legendre.

    We confirm this by computing the sample points and weights over the interval [-1,1] for an odd number of sample points. Since the function operates by then transforming these coordinates and points onto the interval requested, coincidence in the [-1,1] case should ensure coincidence in all other cases. However, we check this for the [0,1] interval with an even number of sample points to be safe
    """
    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_coords, mat_weights = engine.gauss_legendre(a, b, float(n_samples), nargout=2)
    mat_coords = np.array(mat_coords)
    mat_weights = np.array(mat_weights)
    engine.quit()

    py_coords, py_weights = gauss_legendre(n_samples, a, b)

    assert np.all(
        np.isclose(py_coords, mat_coords)
    ), f"Coords different when n_samples={n_samples} over range [{a},{b}]"
    assert np.all(
        np.isclose(py_weights, mat_weights)
    ), f"Weights different when n_samples={n_samples} over range [{a},{b}]"

    return
