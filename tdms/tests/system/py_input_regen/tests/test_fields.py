import os

import numpy as np

from ..fields.detector import gaussian_detector
from ..fields.gauss_pol import gauss_pol_FWHM, gauss_pol_OOES
from ..fields.planewave import planewave_polarised_X, planewave_Z
from ..fields.zero_field import zero_field
from ..interface import Interface
from ..matlab_engine.matlab_engine import create_engine_for_testing
from ..misc_functions.fdtd_bounds import source_grid_coords
from ..timestepping import frequency_vector

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))


def test_planewave_Z() -> None:
    """Test planewave_Z against it's MATLAB counterpart, efield_plane"""
    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)

    # Scalars into meshgrid
    x = y = z = 0.0
    X, Y, Z = np.meshgrid(
        x, y, z, indexing="ij"
    )  # 'ij' indexing is equivalent to matlab's ndgrid function
    # Compute fields both ways
    py_field_scalar = planewave_Z(Z)
    mat_field_scalar = engine.efield_plane(X, Y, Z, nargout=1)
    mat_field_scalar = [np.array(a, dtype=complex) for a in mat_field_scalar]

    assert all(
        np.allclose(py, mat) for py, mat in zip(py_field_scalar, mat_field_scalar)
    )

    # Vectors into meshgrid
    x = y = z = np.linspace(0.0, 1.0, num=5, endpoint=True, dtype=float)
    X, Y, Z = np.meshgrid(x, y, z)

    py_field = planewave_Z(Z)
    mat_field = engine.efield_plane(X, Y, Z, nargout=1)
    mat_field = [np.array(a, dtype=complex) for a in mat_field]

    assert all(np.allclose(py, mat) for py, mat in zip(py_field, mat_field))

    engine.quit()
    return


def test_zero_field() -> None:
    """Test that zero_field returns a list of 3 empty arrays of the expected size."""
    target_shape = (3, 1, 5)

    null_field = zero_field(target_shape)
    print(null_field)

    assert all(component.shape == target_shape for component in null_field)
    assert sum(np.sum(component) for component in null_field) == 0.0

    return


def test_gauss_pol() -> None:
    """
    gauss_pol_OOES is designed to be equivalent to gauss_pol_base with tight = True.
    gauss_pol_FWHM is designed to be equivalent to gauss_pol_base with tight = False.
    """
    # Setup variables
    OOES = 7.5e-6
    FWHM = 21.5e-6
    theta = np.linspace(0.0, 2.0 * np.pi, num=15, endpoint=False)
    phi = np.linspace(0.0, 2.0 * np.pi, num=17, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_OOES = engine.gauss_pol_base(THETA, PHI, True, OOES, nargout=1)
    mat_FWHM = engine.gauss_pol_base(THETA, PHI, False, FWHM, nargout=1)
    engine.quit()

    mat_OOES = np.array(mat_OOES)
    mat_FWHM = np.array(mat_FWHM)

    py_OOES = gauss_pol_OOES(THETA, OOES)
    py_FWHM = gauss_pol_FWHM(THETA, FWHM)

    assert np.all(np.isclose(py_OOES, mat_OOES))
    assert np.all(np.isclose(py_FWHM, mat_FWHM))

    return


def test_gaussian_detector() -> None:
    """gaussian_detector is tested against it's MATLAB counterpart, gaussian_detfun."""
    x = np.linspace(0, 1, num=20, endpoint=True)
    y = np.linspace(1, 2, num=15, endpoint=True)
    X, Y = np.meshgrid(x, y, indexing="ij")
    beta = 25.0 / 36.0
    modal_range = np.arange(0, 31, step=1, dtype=float) * 1.0e-6

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    for m, mode in enumerate(modal_range):
        # MATLAB function takes the index in the allowed modal range,
        # rather than the raw value itself.
        # Also need to offset because of indexing conventions.
        mat_output = engine.gaussian_detfun(X, Y, beta, m + 1)
        mat_output = np.array(mat_output)

        py_output = gaussian_detector(X, Y, beta, mode)

        assert np.all(
            np.isclose(py_output, mat_output)
        ), f"Detector values different at {mode} ({m}-th value in range)"
    engine.quit()

    return


def test_planewave_polarised_X() -> None:
    """planewave_polarised_X should produce the same Exi field as it's MATLAB counterpart, calc_field_tdfield.

    Input file 10 is chosen since arc_10 requires the use of calc_field_tdfield in it's setup.
    """
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_10.m"
    )
    # Where to place the output that MATLAB produces, so we can clean it up later
    DUMP_EIVARS = os.path.abspath(LOCATION_OF_THIS_FILE + "eivars_dump.mat")
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
    refractive_index = 1.35
    dt = 2.0 * 0.95 * delta["x"] * refractive_index / (np.sqrt(2.0) * np.pi * 3e8)
    Nt = 2000
    freq_vector = frequency_vector(Nt, dt)
    f_an = np.arcsin(2.0 * np.pi * 2.997924580105029e8 * dt / (2.0 * 1300e-9)) / (
        np.pi * dt
    )
    interface = Interface(
        i0=[5, 0], i1=[I - 5, 0], j0=[5, 0], j1=[J - 5, 0], k0=[10, 1], k1=[K - 5, 0]
    )
    wavelength_width = 120e-9
    illorigin = np.array(np.floor([I / 2, J / 2, K / 2]))
    sourcemode = "pulsed"

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    engine.calc_field_tdfield(MFILE, DUMP_EIVARS, nargout=0)
    mat_Exi = engine.load(DUMP_EIVARS, "exi")["exi"]
    engine.quit()
    # Cleanup the "output file" that calc_field_fdfield produces with the field values, after we have loaded the values in
    os.remove(DUMP_EIVARS)
    # Convert to numpy arrays
    mat_Exi = np.array(mat_Exi)

    # Call Python planewave_polarised_X, which should produce the same output for Exi
    ex_coords = source_grid_coords(
        delta,
        PML,
        interface,
        illorigin,
        sourcemode,
        I=I,
        J=J,
        K=K,
        requested_components=["Ex"],
    )[0]
    dz = delta["z"] / 2

    py_Exi = planewave_polarised_X(
        ex_coords, freq_vector, refractive_index, f_an, wavelength_width, dz=dz
    )

    # py_Exi returns an array that retains the shape of the coordinate array (IE no flattening takes place)
    assert np.all(
        py_Exi.shape
        == (ex_coords["x"].size, ex_coords["y"].size, ex_coords["z"].size, Nt)
    )
    # Assert that, after flattening dimensions of size 1, the arrays have the same shape.
    # This needs to be done because the matlab function assumes that the input lies exclusively on the K0 plane, and thus collapses dimensions automatically.
    assert np.squeeze(py_Exi).shape == np.squeeze(mat_Exi).shape
    # Assert that the contents of the arrays are the same
    py_Exi = np.squeeze(py_Exi)
    mat_Exi = np.squeeze(mat_Exi)
    assert np.all(np.isclose(py_Exi, mat_Exi))

    return
