from typing import Literal

import numpy as np


def multi_layer(
    zvec: np.ndarray,
    nvec: np.ndarray,
    wavelength_in_air: float,
    theta0,
    polarisation: Literal["TE", "TM"] = "TE",
) -> np.ndarray:
    """
    Compute the characteristic matrix of a multilayer structure.

    The input array nvec should have one more element than the zvec input.

    :param zvec: Vector defining the z-locations of the layers of the (unrotated) stratified medium.
    :param nvec: Refractive indices of the layers.
    :param wavelength_in_air: Wavelength of light in air.
    :param theta0: Angle of incidence of the plane wave, in radians.
    :param polarisation: The polarisation of the source. 'TE' means that the electric vector is perpendicular to the plane in which the wave propagates. 'TM' is used when the electric field is parallel to, or in, the plane in which the wave propagates.
    :returns: The characteristic matrix of the multiplayer structure.
    """
    # Check dimensionality of input
    if nvec.size != (zvec.size + 1):
        raise RuntimeError(
            f"nvec should have one more element than zvec, but they have {nvec.size} and {zvec.size} elements respectively."
        )

    # Constant values
    theta_vec = np.arcsin(nvec[0] * np.sin(theta0) / nvec)
    cos_theta = np.cos(theta_vec)
    # Speedup to save recomputing multiple times
    intermediate_product = nvec * cos_theta

    if nvec.size > 2:
        beta = (
            intermediate_product[1:-1] * np.diff(zvec) * 2.0 * np.pi / wavelength_in_air
        )

    if polarisation == "TE":
        denominator = intermediate_product[:-1] + intermediate_product[1:]
        t_fresnel = 2.0 * intermediate_product[:-1] / denominator
        r_fresnel = -np.diff(intermediate_product) / denominator
    else:  # polarisation == 'TM'
        denominator = nvec[1:] * cos_theta[:-1] + nvec[:-1] * cos_theta[1:]
        t_fresnel = 2.0 * intermediate_product[:-1] / denominator
        r_fresnel = (
            nvec[1:] * cos_theta[:-1] - nvec[:-1] * cos_theta[1:]
        ) / denominator

    # Define characteristic/scattering matrix
    S = 1 / t_fresnel[0] * np.array([[1.0, r_fresnel[0]], [r_fresnel[0], 1.0]])
    for i in range(1, nvec.size - 1):
        # Note that beta is not defined if nvec has 2 or fewer elements.
        # However in this case, this loop doesn't execute.
        P = np.diag(np.exp([-1.0j * beta[i - 1], 1.0j * beta[i - 1]]))
        L = 1 / t_fresnel[i] * np.array([[1, r_fresnel[i]], [r_fresnel[i], 1]])
        # S = S * P * L
        S = np.matmul(np.matmul(S, P), L)

    return S
