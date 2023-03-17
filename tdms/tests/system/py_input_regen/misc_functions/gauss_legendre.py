from typing import Tuple

import numpy as np

NEWTON_THRESHOLD = np.spacing(1.0)


def gauss_legendre(
    n_sample_points: int, a: float = -1, b: float = 1
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the absciassa (sample coordinates) and weights for performaing numerical integration within the interval (a,b), using n_sample_points via Gauss-Legendre quadrature.

    From these points, the integral
        F = \int_a^b f(x) dx,
    can be evaluated as
        F = sum( f(sample_coords) .* sample_weights ).

    :param n_sample_points: The number of sample points to use.
    :param a,b: Endpoints of the interval of integration.
    :returns: Tuple[sample_coords, sample_weights], as detailed above.
    """
    # There will be n_sample_points roots (coordinates) and weights, which are symmetric
    sample_coords = np.zeros((n_sample_points,), dtype=float)
    sample_weights = np.zeros((n_sample_points,), dtype=float)

    # Compute sample_coords and sample_weights for the interval (-1,1)
    for i in range(0, int(np.ceil(n_sample_points / 2))):
        # Tracker for Newton convergence
        converged = False

        # Use an asymptotic expansion ("Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables, Milton Abramowitz and Irene A. Stegun, Dover, New York, 1964, Eq 8.10.8) to estimate the value of the root
        x0 = np.cos((i + 3 / 4) * np.pi / (n_sample_points + 1 / 2))

        while not converged:
            # Now we need to evaluate the value of the Legendre polynomial and its derivative, in order to form a Taylor series expansion of order 1 at x0.

            # Evaluate the Legendre polynomial at x0 (Abramowitz and Stegun, Eq 8.5.3)
            p0 = 1.0
            p1 = x0
            for j in range(2, n_sample_points + 1):
                p2 = ((2 * j - 1) * x0 * p1 - (j - 1) * p0) / j
                p0 = p1
                p1 = p2
            # Evaluate the derivative (Abramowitz and Stegun, Eq 8.5.4)
            p2_prime = n_sample_points * (x0 * p1 - p0) / (x0**2 - 1)

            # Determine if convergence has been achieved
            xml = x0
            x0 -= p1 / p2_prime
            converged = not (np.abs(x0 - xml) > NEWTON_THRESHOLD)

        # Having obtained this sample point, compute the weight and use symmetry to compute the correspoding weight and sample point
        # Coords have antisymmetry
        sample_coords[i] = -x0
        sample_coords[-(i + 1)] = x0
        # Weights have symmetry
        sample_weights[i] = 2.0 / ((1.0 - x0**2) * p2_prime**2)
        sample_weights[-(i + 1)] = sample_weights[i]

    # Transform the coordinates and weights from the interval (0,1) to (a,b), via the transformation
    # x = e * tau + f,
    # mapping tau = -1 to x = a, and tau = 1 to x = b.
    # This results in
    # e = (b - a) / 2,
    # f = (b + a) / 2
    e = (b - a) / 2
    f = (b + a) / 2
    sample_coords = e * sample_coords + f
    sample_weights = e * sample_weights

    return sample_coords, sample_weights
