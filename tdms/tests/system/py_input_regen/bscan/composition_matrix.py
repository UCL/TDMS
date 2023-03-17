from typing import Literal

import numpy as np


def build_composition_matrix(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    obstacle_shape: Literal["fs", "cyl", "sph", "sc"],
    obstacle_radius: float = 1.0,
) -> np.ndarray:
    """
    Create the composition matrix for the spatial grid provided and obstacle specified.

    The composition matrix is a 2D array of floats (COULD BE INTS W/ TDMS API CHANGES) whose rows are of the format [i, j, k, 1] where (i,j,k) is the index of a Yee cell whose centre is contained within the spatial object.

    :param x,y,z: (N,) numpy arrays whose Cartesian product are the coordinates of the Yee cell centres.
    :param obstacle_shape: The shape of the scattering obstacle, either 'fs' (free-space), 'cyl' (cylindrical), 'sph' (spherical), or 'sc' (point-source).
    :param obstacle_radius: The radius of the scattering obstacle, ignored by free-space or point-source obstacles.
    :returns: The composition matrix, detailed above.
    """
    # Free space composition matrix is just the empty array
    if obstacle_shape == "fs":
        return np.empty((0, 0), dtype=float)

    # If defining a cylindrical object, we need to drop the y-dimension
    if obstacle_shape == "cyl":
        yy = 0
    else:
        yy = y
    X, Y, Z = np.meshgrid(x, yy, z, indexing="ij")

    I = np.zeros_like(X, dtype=float)
    # Set all Yee cells within the obstacle to have an index of 1
    if obstacle_shape == "cyl":
        I[X**2 + Z**2 < obstacle_radius**2] = 1
        I[-3:, 0, :] = 0
        I[:, 0, -3:] = 0
    elif obstacle_shape == "sph":
        I[X**2 + Y**2 + Z**2 < obstacle_radius**2] = 1
        I[-3:, 0, :] = 0
        I[:, 0, -3:] = 0
    elif obstacle_shape == "sc":
        obstacle_cells = np.logical_and(np.logical_and(X == 0, Y == 0), Z == 0)
        I[obstacle_cells] = 1
    else:
        raise ValueError(f"Error: {obstacle_shape} is not a valid obstacle shape.")

    # Format composition matrix as a 2D array whose rows are [i, j, k, 1]; where i,j,k is the index of a Yee cell within the scattering obstacle
    indices = np.argwhere(I)
    composition_matrix = np.c_[indices, np.ones((indices.shape[0],), dtype=float)]

    return composition_matrix
