from copy import deepcopy
from typing import Literal, Tuple

import numpy as np

from ..interface import Interface
from .yee_cell_position import yee_position


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


def source_grid_coords(
    delta: dict[str, float],
    PML: dict[str, int],
    interface: Interface,
    illorigin: np.ndarray,
    sourcemode: Literal["pulsed", "steadystate"],
    I: int,
    J: int,
    K: int,
    requested_plane: Literal["I0", "I1", "J0", "J1", "K0", "K1"] = "K0",
    requested_components: list[str] = ["Ex", "Ey"],
    z_launch=0.0,
) -> Tuple[dict[str, np.ndarray]]:
    """Computes the spatial coordinates of the source field, on the requested plane and for the requested field components.

    Outputs are returned as a tuple of dicts. Each dict corresponds to one field component that has been provided in requested_components.
    Dicts are returned in alphabetical order. A given dictionary then has three keys; 'x', 'y', 'z', which are arrays containing the spatial coordinates of that field component on the plane.

    Providing plane = 'K0' and requested_components = ['Hx', 'Ez'] for example, returns a tuple of two dictionaries. These dictionaries correspond to the 'Ez' and 'Hx' components, in that order. The values of the first dict are the spatial coordinates of the Ez component, and the second dict those of the Hx component.

    We use a total/scattered field formulation. For k < K1 we have scattered, and for k >= K1 we have total. The field is coupled into the system by means of update equations to maintain consistency.

    :param delta: Yee cell extents in each of the coordinate dimensions.
    :param PML: The number of Yee cells in the Perfectly Matched Layer on each face of the spatial grid.
    :param interface: Specifications for the interface conditions on each of the planes {I,J,K}{0,1}.
    :param illorigin: Coordinate corresponding to the origin of the illumination/source field.
    :param sourcemode: Whether the source field is introduced at steady state, or pulsed.
    :param I,J,K: Number of Yee cells in the x,y,z directions respectively.
    :param requested_plane: The plane or interface on which we want the coordinates of the field components.
    :param requested_components: The field components that we want to obtain spatial coordinates for.
    :param z_launch: Additional "height" above the interface of the source field.
    :returns: Tuple of dicts. The order of the dicts corresponds to the requested field components, sorted into alphabetical order. Each dictionary has fields 'x', 'y', and 'z' containing the spatial coordinates on the plane of the given component.
    """
    ## Take deepcopies of objects we may need to adjust
    interface = deepcopy(interface)
    illorigin = deepcopy(illorigin)
    # Source origin always needs to be cast to global coordinates
    illorigin += np.array([PML["Dxl"], PML["Dyl"], PML["Dzl"]])

    ## Validate the interface_conditions that have been provided
    interface.validate_against_n_cells(I, J, K)

    ## Prepare the interface, and change indexing system if necessary
    if sourcemode == "pulsed":
        # We only use the K0 plane in this case, set defaults for the other components
        interface["I0"][0] = interface["J0"][0] = 1
        interface["I1"][0] = I + PML["Dxl"] + PML["Dxu"] + 1
        interface["J1"][0] = J + PML["Dyl"] + PML["Dyu"] + 1
        for plane in ["I0", "I1", "J0", "J1", "K1"]:
            interface[plane][1] = 0
        if (not interface.non_zero_interface("K0")) and K != 0:
            raise RuntimeError(
                f"Error: Running in pulsed mode with no K0 interface condition - no point running simulation."
            )
    elif sourcemode == "steadystate":
        # We need to use global coordinates
        interface.change_coords("global", PML)
    else:
        # This isn't a valid source mode!
        raise ValueError(f"Error: {sourcemode} is not a valid source mode.")

    ## Extract the field coordinates on the requested plane
    # Infer the cell indices that we need to lookup spatial coordinates for.
    i_source = interface.cell_range("I") - illorigin[0]
    j_source = interface.cell_range("J") - illorigin[1]
    k_source = interface.cell_range("K") - illorigin[2]
    if "I" in requested_plane:
        i_source = interface[requested_plane][0] - illorigin[0]
    elif "J" in requested_plane:
        j_source = interface[requested_plane][0] - illorigin[1]
    elif "K" in requested_plane:
        k_source = interface[requested_plane][0] - illorigin[2]
    # The Cartesian product i_source \otimes j_source \otimes k_source now provides the set of (i,j,k) Yee cell indices we need to lookup.
    # Note that one of i_source, j_source, k_source is a scalar, corresponding to the dimension the plane lies perpendicular to.

    # Now fetch the coordinates for each of the components that were requested.
    coords = dict()
    for component in requested_components:
        coords[component] = dict()
        if interface.non_zero_interface(requested_plane):
            # Fetch the coordinates of this field component on the requested plane
            (
                coords[component]["x"],
                coords[component]["y"],
                coords[component]["z"],
            ) = yee_position(i_source, j_source, k_source, delta, component)
            # Adjust z coords for possible height offset
            coords[component]["z"] += z_launch
        else:
            # This component's coordinates were not requested, just return empty arrays for the components
            coords[component]["x"] = coords[component]["y"] = coords[component][
                "z"
            ] = np.array([], dtype=float)

    return tuple(coords.values())
