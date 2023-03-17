from typing import Tuple

import numpy as np

from .gaussians import gaussian_1D


def gauss_pol_OOES(
    theta: np.ndarray,
    OOES: float = 5.0e-6,
    wavelength: float = 1300.0e-9,
    relative_refindex: float = 1.35,
) -> Tuple[np.ndarray, np.ndarray]:
    """Creates a Gaussian pulse G(theta) of the form
        G(theta) = exp( - (sin(theta) * W)^2 ),
    where the constant
        W = wavenumber * relative_refindex * OOES / 4.

    :param theta: Angular values to evaluate at
    :param OOES: FIXME
    :param wavelength: Wavelength of light in free space
    :param relative_refindex: Relative refractive index of the material
    :returns: Tuple, (G(theta), G(theta))
    """
    wavenumber = 2.0 * np.pi / wavelength
    W = wavenumber * relative_refindex * OOES / 4.0

    wx = gaussian_1D(np.sin(theta), variance=1.0 / W)

    return wx, wx


def gauss_pol_FWHM(
    theta: np.ndarray,
    FWHM: float = 25.0e-6,
    wavelength: float = 1300.0e-9,
    relative_refindex: float = 1.35,
) -> Tuple[np.ndarray, np.ndarray]:
    """Creates a Gaussian pulse G(theta) of the form
        G(theta) = exp( - (sin(theta) * W)^2 ),
    where the constant
        W = wavenumber * relative_refindex * FWHM / (2. * sqrt(2 * ln(2)))

    :param theta: Angular values to evaluate at
    :param FWHM: FIXME
    :param wavelength: Wavelength of light in free space
    :param relative_refindex: Relative refractive index of the material
    :returns: Tuple, (G(theta), G(theta))
    """
    wavenumber = 2.0 * np.pi / wavelength
    W = FWHM * wavenumber * relative_refindex / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    wx = gaussian_1D(np.sin(theta), variance=1.0 / W)

    return wx, wx
