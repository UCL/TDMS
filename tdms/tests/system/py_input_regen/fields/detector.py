import numpy as np

from .gaussians import gaussian_2D


def gaussian_detector(
    X: np.ndarray, Y: np.ndarray, beta: float, modal_value: float
) -> np.ndarray:
    """
    Defines a 2D-Gaussian detector function with x_mean modal_value and variance equal to 9.2e-6 / beta.

    See also: gaussians.gaussian_2D
    """
    MFD = 9.2e-6
    variance = MFD / beta

    return gaussian_2D(
        X, Y, x_mean=modal_value, x_variance=variance, y_variance=variance
    )
