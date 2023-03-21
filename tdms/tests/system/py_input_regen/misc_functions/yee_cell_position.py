from typing import Literal, Tuple, Union

import numpy as np


def yee_position(
    i: Union[int, np.ndarray],
    j: Union[int, np.ndarray],
    k: Union[int, np.ndarray],
    delta: dict[str, float],
    component: Literal["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"] = "Ex",
) -> Tuple[
    Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray]
]:
    """Calculate the position in Cartesian space of the field component associated to the Yee cell with index (i,j,k).

    E-field components are offset by 0.5 * the Yee cell extent in that dimension.
    H-field components are offset by 0.5 * the Yee cell extent in the OTHER two dimensions.

    :param i,j,k: Triples that form (i,j,k) indices of Yee cells.
    :param delta: Dict with keys 'x','y','z' that correspond to the Yee cell extent in the respective dimensions.
    :param component: One of 'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'. The field component we are examining.
    """
    # Components align with the Yee cell centre, save for the field component we are examining
    x = i * delta["x"]
    y = j * delta["y"]
    z = k * delta["z"]

    # Add offsets depending on the component we are examining
    if component == "Ex":
        x += 0.5 * delta["x"]
    elif component == "Ey":
        y += 0.5 * delta["y"]
    elif component == "Ez":
        z += 0.5 * delta["z"]
    elif component == "Hx":
        y += 0.5 * delta["y"]
        z += 0.5 * delta["z"]
    elif component == "Hy":
        x += 0.5 * delta["x"]
        z += 0.5 * delta["z"]
    else:
        # component == 'Hz'
        x += 0.5 * delta["x"]
        y += 0.5 * delta["y"]

    return x, y, z
