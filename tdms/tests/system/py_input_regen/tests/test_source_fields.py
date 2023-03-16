import os

import numpy as np

from ..matlab_engine.matlab_engine import create_engine_for_testing
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
