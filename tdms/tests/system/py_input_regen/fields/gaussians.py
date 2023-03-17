import numpy as np


def gaussian_1D(X: np.ndarray, mean: float = 0.0, variance: float = 1.0) -> np.ndarray:
    """
    Evaluates the scalar Gaussian G(x),
        G(x) = exp( -((x-x_mean)/x_variance)**2 )
    over the coordinates provided.

    :param X: Coordinates to evaluate Gaussian at.
    :param mean: Coordinate translations.
    :param variance: Scale factor.
    """
    XX = ((X - mean) / variance) ** 2

    return np.exp(-XX)


def gaussian_2D(
    X: np.ndarray,
    Y: np.ndarray,
    x_mean: float = 0.0,
    x_variance: float = 1.0,
    y_mean: float = 0.0,
    y_variance: float = 1.0,
) -> np.ndarray:
    """
    Evaluates the scalar Gaussian G(x,y),
        G(x,y) = exp( - ( ((x-x_mean)/x_variance)**2 + ((y-y_mean)/y_variance)**2 ) )
    over the coordinates provided.

    :param X,Y: Coordinates to evaluate Gaussian at.
    :param x_mean,y_mean: Shifts to the x and y coordinate values, respectively.
    :param x_variance,y_variance: Scalar factors for the x and y coordinates, respectively.
    """
    XX = ((X - x_mean) / x_variance) ** 2
    YY = ((Y - y_mean) / y_variance) ** 2

    return np.exp(-(XX + YY))
