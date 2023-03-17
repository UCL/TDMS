import os

import numpy as np

from ..matlab_engine.matlab_engine import create_engine_for_testing
from ..source_fields.gauss_pol import gauss_pol_FWHM, gauss_pol_OOES
from ..source_fields.planewave import planewave_Z
from ..source_fields.zero_field import zero_field

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))


def test_planewave_Z() -> None:
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
