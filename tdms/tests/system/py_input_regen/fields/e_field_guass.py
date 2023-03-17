import numpy as np


def e_field_guass(
    X: np.ndarray, Y: np.ndarray, Z: np.ndarray, guass_pol_method, **args
):
    """ """

    # Setup the output field
    E = [np.zeros_like(X), np.zeros_like(X)]

    # Constants that really should be input arguments
    wavelength = 1300e-9
    relative_refindex = 1.35
