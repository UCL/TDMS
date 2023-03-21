import os
from typing import Literal

import numpy as np
import pytest
from pytest_check import check

from ..interface import Interface
from ..matlab_engine.matlab_engine import create_engine_for_testing
from ..misc_functions.fdtd_bounds import fdtd_grid_coords, source_grid_coords
from ..misc_functions.gauss_legendre import gauss_legendre
from ..misc_functions.multi_layer import multi_layer
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


@pytest.mark.parametrize("polarisation", ["TM", "TE"])
def test_multi_layer(polarisation: Literal["TE", "TM"]) -> None:
    """multi_layer should be equivalent to it's MATLAB counterpart, also named multi_layer.

    We test (via parameterisation) both the TE and TM polarisation settings.
    """
    wavelength = 1300e-9
    # Refractive indices to trial
    nvec = np.array([1.0, 1.35, 2.25])
    # Angles to sweep over in one test
    incidence_angles = np.linspace(0, 2 * np.pi, num=10, endpoint=False)
    # Some arbitrary "layer heights" to use
    zvec = np.arange(nvec.size - 1, dtype=float) * 1.5

    # Assert that nvec must have one more element that zvec
    with pytest.raises(RuntimeError):
        multi_layer(np.arange(nvec.size), nvec, wavelength, 0.0, polarisation)

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    with check:
        for theta0 in incidence_angles:
            mat_matrix = engine.multi_layer(
                zvec, nvec, wavelength, theta0, polarisation, nargout=1
            )
            mat_matrix = np.array(mat_matrix)

            py_matrix = multi_layer(zvec, nvec, wavelength, theta0, polarisation)

            assert np.all(
                np.isclose(py_matrix, mat_matrix)
            ), f"Composition matrices differ at angle {theta0}"
    engine.quit()

    return


def test_fdtd_bounds() -> None:
    """fdtd_grid_coords should be equivalent to the MATLAB method fdtd_bounds when the inputs from the input file match those provided to the Python function.

    MATLAB function takes an input file as argument. Choice is somewhat arbitrary, so we use input_file_08.m
    """
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_08.m"
    )
    # Variables in the input_file that Python will need to read
    delta = {
        "x": 1300.0e-9 / 8.0,
        "y": 1300.0e-9 / 8.0,
        "z": 1300.0e-9 / 8.0,
    }
    I = 32
    J = 32
    K = 16
    illorigin = np.array(np.floor([I / 2, J / 2, K / 2]))
    # z_launch is 0, no test exists where this is non-zero

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    # MATLAB function also reads the wavelength from the input file, which we don't need here
    x_mat, y_mat, z_mat, _ = engine.fdtd_bounds(MFILE, nargout=4)
    x_mat = np.array(x_mat)
    y_mat = np.array(y_mat)
    z_mat = np.array(z_mat)
    engine.quit()

    x_py, y_py, z_py = fdtd_grid_coords(delta, I, J, K, illorigin)

    with check:
        assert np.all(np.isclose(x_py, x_mat)), f"Grid coordinates in x differ."
        assert np.all(np.isclose(y_py, y_mat)), f"Grid coordinates in y differ."
        assert np.all(np.isclose(z_py, z_mat)), f"Grid coordinates in z differ."

    return


def test_source_grid_coords() -> None:
    """
    source_grid_coords should produce the same coordinates for the Ex and Ey components on the K0 plane as the MATLAB equivalent getsourcecoords.

    We test against input_file_10.m because this test requires calc_tdfield, thus will require this method to be called.
    """
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_10.m"
    )
    # Variables in the input_file that Python will need to read
    delta = {
        "x": 1300.0e-9 / 6.0,
        "y": 1300.0e-9 / 6.0,
        "z": 1300.0e-9 / 6.0,
    }
    PML = {
        "Dxl": 10,
        "Dxu": 10,
        "Dyl": 0,
        "Dyu": 0,
        "Dzl": 10,
        "Dzu": 10,
    }
    I = 64
    J = 0
    K = 64
    interface = Interface(
        i0=[5, 0], i1=[I - 5, 0], j0=[5, 0], j1=[J - 5, 0], k0=[10, 1], k1=[K - 5, 0]
    )
    illorigin = np.array(np.floor([I / 2, J / 2, K / 2]))
    sourcemode = "pulsed"
    # z_launch is 0, no test exists where this is non-zero

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    # We only produce the first two outputs, ex_coords and ey_coords, with our converted source_grid_coords function.
    mat_ex, mat_ey = engine.getsourcecoords(MFILE, nargout=2)
    engine.quit()

    # By default, source_grid_coords extracts coords on K0, for the Ex and Ey coordinates
    py_ex, py_ey = source_grid_coords(
        delta, PML, interface, illorigin, sourcemode, I=I, J=J, K=K
    )

    for axis in ["x", "y", "z"]:
        assert np.all(
            np.isclose(py_ex[axis], mat_ex[axis])
        ), f"Ex coords along axis {axis} do not match."
        assert np.all(
            np.isclose(py_ey[axis], mat_ey[axis])
        ), f"Ey coords along axis {axis} do not match."

    return
