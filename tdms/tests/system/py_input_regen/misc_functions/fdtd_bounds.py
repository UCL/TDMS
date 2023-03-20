from typing import Tuple

import numpy as np


def fdtd_grid_coords(
    delta: dict[str, float],
    I: int,
    J: int,
    K: int,
    illorigin: np.ndarray,
    z_launch: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute 1D-vectors x, y, z which contain the coordinates of the FDTD grid in the respective axial direction. That is, the Cartesian product x \otimes y \otimes z is the set of all gridpoints (x_i, y_j, z_k) in the FDTD grid.

    :param delta: Yee cell extents in each of the coordinate dimensions.
    :param I,J,K: Number of Yee cells in the x,y,z directions respectively.
    :param illorigin: Coordinate corresponding to the origin of the illumination/source field.
    :param z_launch: Additional "height" above the interface of the source field.
    :returns: (x, y, z), vectors containing the coordinates of the FDTD grid.
    """
    x = (np.arange(1, I + 1, dtype=float) - illorigin[0]) * delta["x"]
    y = (np.arange(1, J + 1, dtype=float) - illorigin[1]) * delta["y"]
    z = (np.arange(1, K + 1, dtype=float) - illorigin[2]) * delta["z"]
    # Adjust for source field "elevation"
    z += z_launch

    return x, y, z
