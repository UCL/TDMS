from typing import Union

import numpy as np
from numpy.fft import fft

from ..physical_constants import C
from .gaussians import gaussian_1D, half_width_half_maximum


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


def planewave_polarised_X(
    coords: dict[str, np.ndarray],
    freq_vector: np.ndarray,
    refractive_index: float,
    f_an: float,
    wavelength_width: float,
    dz: float = 0.0,
) -> np.ndarray:
    """Defines a planewave linearly polarised in the x-direction, and evaluates it over the spatial grid provided. The range of wavelengths to evaluate over is determined by frequency_vector(Nt, dt).

    The functional form of the field component in real space is
        E_{i}(x, y, z, lambda)
            = -2i\pi HWHM G(C/lambda) * exp((2i\pi to_l C/lambda),
    where
        G(theta) = exp( - \pi HWHM^2 (theta - f_an)^2 ),
        to_l = HWHM * sqrt(log(1e8)/\pi),
    and HWHM is the half-width half-maximum of a Gaussian pulse with variance
        sigma = 2 * ref_index * lambda_an^2 / (C * wavelength_width * sqrt(pi))
            = 2 * C / (wavelength_width * ref_index * f_an^2 * sqrt(pi))
    This field then has a Fourier transform taken along the lambda (wavelength) axis, and the real part of the resulting Fourier-space field is returned as the (time-domain) field.

    The output is a floating point numpy array of the field evaluations over the spatial grid and wavelength range provided. It has shape (coords['x'].size, coords['y'].size, coords['z'].size, Nt)

    :param coords: The spatial coordinates over which to evaluate the field. Should be a dict with keys 'x', 'y', 'z', that have 1D arrays as their values, such that the Cartesian product of these arrays forms the set of spatial coordinates to evaluate the field over.
    :param freq_vector: The frequencies corresponding to the wavelengths the field should be evaluated over; freq_vector = c / lambda.
    :param refractive_index: Refractive index of the medium.
    :param f_an: FIXME (average?)
    :param wavelength_width: FIXME
    :param dz: Offset to apply to the z-coordinates provided.
    :returns: E_{i}(x, y, z, lambda) evaluated across the spatial grid and wavelength range.
    """
    x, y, z = coords["x"], coords["y"], coords["z"] - dz
    # Setup Gaussian-related constants
    sigma = 2 * C / (wavelength_width * refractive_index * f_an**2 * np.sqrt(np.pi))
    HWHM = half_width_half_maximum(sigma)
    to_l = HWHM * np.sqrt(np.log(1.0e8) / np.pi)

    # Compute G(C/lambda)
    X, Y, Z, C_OVER_LAMBDA = np.meshgrid(x, y, z, freq_vector, indexing="ij")
    G = gaussian_1D(C_OVER_LAMBDA, mean=f_an, variance=1.0 / (HWHM * np.sqrt(np.pi)))

    # Compute the pre-FFT field
    prefactor = -2.0j * HWHM
    phase_factor = np.exp(2.0j * np.pi * to_l * C_OVER_LAMBDA)
    field_values = prefactor * G * phase_factor

    # Fourier transform along the lambda-axis
    field_values = fft(field_values, axis=3)
    # Rescale by the frequency difference
    field_values *= freq_vector[1] - freq_vector[0]
    # Take real part
    field_values = np.real(field_values)
    # Shift first fft component to end of each (hyper-) column along the wavelength-axis
    field_values = np.roll(field_values, -1, axis=3)

    return field_values
