from typing import Tuple

import numpy as np


def gauss_pol_OOES(
    th: np.ndarray,
    OOES: float = 5.0e-6,
    wavelength: float = 1300.0e-9,
    relative_refindex: float = 1.35,
) -> Tuple[np.ndarray, np.ndarray]:
    wavenumber = 2.0 * np.pi / wavelength

    W = 4.0 / (wavenumber * relative_refindex * OOES)
    wx = np.exp(-((np.sin(th) / W) ** 2))

    return wx, wx


def gauss_pol_FWHM(
    th: np.ndarray,
    FWHM: float = 25.0e-6,
    wavelength: float = 1300.0e-9,
    relative_refindex: float = 1.35,
) -> Tuple[np.ndarray, np.ndarray]:
    wavenumber = 2.0 * np.pi / wavelength

    W = FWHM * wavenumber * relative_refindex / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    wx = np.exp(-((W * np.sin(th)) ** 2))

    return wx, wx
