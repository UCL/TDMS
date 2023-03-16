from typing import Union

import numpy as np


def planewave_Z(
    Z: np.ndarray,
    wavelength: float = 1300.0e-9,
    relative_refindex: float = 1.35,
    dz: Union[float, None] = None,
) -> list[np.ndarray]:
    """Defines a plane-wave propagating in the Z-direction, over a the spatial grid provided.

    :param Z: z-coordinates of the spatial grid the field is passed over
    :param wavelength: Wavelength of light in the medium
    :param relative_refindex: Relative refractive index of the medium
    :param dz: Offset to apply to z-coordinates to account for illumination source
    :returns: The (complex) plane-wave field propagating in the z-direction
    """
    # Set Z-component offset to default if dz was not passed in, otherwise use default
    if dz is None:
        delta_z = wavelength / 10
    else:
        delta_z = dz

    wavenumber = 2.0 * np.pi / wavelength
    # Compute the field
    field = [np.zeros_like(Z, dtype=complex), np.zeros_like(Z, dtype=complex)]
    field[0] = 2.0 * np.exp(1.0j * wavenumber * relative_refindex * (Z - delta_z / 2.0))

    return field
