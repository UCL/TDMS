from collections import UserDict
from typing import Any, Literal, Union

import numpy as np

# Keys that the Interface object will possess
INTERFACE_PLANES: list[str] = ["I0", "I1", "J0", "J1", "K0", "K1"]


def check_valid_key(key_to_be: Any) -> None:
    """Check that a valid key has been passed in."""
    if key_to_be not in INTERFACE_PLANES:
        raise KeyError(f"{key_to_be} cannot be assigned to an Interface object")
    return


def check_valid_item(item_to_be: Any) -> np.ndarray:
    """Check that the object passed can be converted to an interface member. If it can, return the object as cast to the expected type."""
    # First, try casting the item
    cast_item = np.array(item_to_be, dtype=int)

    # Check that we have two elements and the expected shape
    if cast_item.size != 2:
        raise ValueError(
            f"Interface members have 2 elements, but recieved {cast_item.size}"
        )
    elif cast_item.ndim != 1:
        raise ValueError(
            f"Interface members have one dimension, but recieved and object with {cast_item.ndim}"
        )
    else:
        cast_item.resize((2,))

    return cast_item


class Interface(UserDict):
    """Defines the boundary (interface) conditions at each of the 6 faces of the FDTD grid, I0, I1, J0, J1, K0, K1. Planes can be referenced with lower-case letters for convenience.

    Planes that are not provided on initialisation are set to impose no boundary conditions, but can be manually updated later. Planes are stored as integer np.ndarrays of shape (2,).

    The first element of each entry corresponds to the Yee cell index that the plane lies in, so for example the I0 plane having 0th entry 2 means that the interface I0 consists of Yee cells with indices (2, j, k). It should be at least 1.

    If both the 0 and 1 plane possess interface conditions, the cell index of the 0th plane should be no greater than that of the 1st plane.

    CAN BE SIMPLIFED - ONLY INCLUDE KEYS THAT HAVE NON-ZERO INTERFACES, SINCE THE [1]ST ELEMENTS ARE JUST USED AS FLAGS. SAVE TIME MAYBE?
    """

    # Current Yee cell indexing system this object is using
    # Global indexing starts from the PML cells, local indexing does not include these cells.
    in_coords: Literal["global", "local"]

    def __init__(self, in_coords: Literal["global", "local"] = "local", **kwargs):
        super().__init__(kwargs)

        # Co-ordinate system we initialised in
        self.in_coords = in_coords

        # Add zero-interface members for functionality
        for plane in INTERFACE_PLANES:
            if plane not in self.keys():
                self[plane] = np.array([0, 0], dtype=int)
        return

    def __setitem__(self, key: str, item: np.ndarray) -> None:
        check_valid_key(key.upper())
        cast_item_safely = check_valid_item(item)
        return super().__setitem__(key.upper(), cast_item_safely)

    def __getitem__(self, key: str) -> np.ndarray:
        return super().__getitem__(key.upper())

    def __str__(self) -> str:
        print_msg = f"Interface object in {self.in_coords} indexing convention, with non-zero conditions:"
        for plane, i_condition in self.items():
            if np.all(i_condition):
                print_msg += f"\n\t{plane}: [{i_condition[0]}, {i_condition[1]}]"
        return print_msg

    def change_coords(
        self, new_coords: Literal["global", "local"], PML: dict[str]
    ) -> None:
        """Convert the cell indices provided in the boundary conditions to the requested indexing method.

        Global indexing includes the PML cells, and begins counting from the first PML cell in each axial direction.
        Local indexing excludes the PML cells.

        :param new_coords: The coordinate system to transform into.
        :param PML: The extent of the PML on each face of the spatial grid.
        """
        if self.in_coords != new_coords:
            # Change coordinates appropriately
            # Global to local: all interfaces - PML (D*l)
            # Local to global: all interfaces + PML (D*l)
            if new_coords == "global":
                sign = 1
            elif new_coords == "local":
                sign = -1
            else:
                raise RuntimeError(f"{new_coords} is not a valid indexing system!")

            self["I0"] += sign * PML["Dxl"]
            self["I1"] += sign * PML["Dxl"]
            self["J0"] += sign * PML["Dyl"]
            self["J1"] += sign * PML["Dyl"]
            self["K0"] += sign * PML["Dzl"]
            self["K1"] += sign * PML["Dzl"]
            # Update variable tracking which coordinates we're currently in
            self.in_coords = new_coords

        return

    def cell_range(self, plane: Literal["I", "J", "K"]) -> np.ndarray:
        """Return the Yee cell indices between the start and end of the interface.

        EG: If plane = 'J', we return the integers between self['J0'][0] and self['J1'][0] INCLUSIVE.

        :param plane: The plane of the interface, either I, J, or K.
        :returns: The Yee cell indices between the edges of the plane.
        """
        P0 = plane + "0"
        P1 = plane + "1"
        if not (self.non_zero_interface(P0) and self.non_zero_interface(P1)):
            raise ValueError(
                f"One of these planes has no boundary condition! {P0}: {self.non_zero_interface(P0)}, {P1}: {self.non_zero_interface(P1)}"
            )
        return np.arange(self[P0][0], self[P1][0] + 1, step=1, dtype=int)

    def non_zero_interface(
        self, plane: Literal["I0", "I1", "J0", "J1", "K0", "K1"]
    ) -> bool:
        """Determine whether the interface condition at the plane provided is non-zero.

        :param plane: The plane of interest
        :returns: True if the plane has a non-zero interface condition. Zero otherwise.
        """
        return bool(self[plane][1])

    def validate_against_n_cells(
        self,
        I: Union[int, None] = None,
        J: Union[int, None] = None,
        K: Union[int, None] = None,
    ):
        """Validate the interface conditions, given the number of cells in each direction of the spatial grid.

        :param I,J,K: Number of Yee cells in the x,y,z directions respectively. Values set to None are skipped.
        """
        planes_to_cells = {"I": I, "J": J, "K": K}
        for direction, n_cells in planes_to_cells.items():
            P0 = direction + "0"
            P1 = direction + "1"
            if (
                (n_cells is not None)
                and self.non_zero_interface(P0)
                and self.non_zero_interface(P1)
            ):
                if self[P0][0] < 1:
                    raise ValueError(f"{P0} interface out of bounds: {self[P0][0]} < 1")
                elif self[P0][0] > self[P1][0]:
                    raise ValueError(
                        f"Interface {P1} has lower cell index ({self[P1][0]}) than {P0} ({self[P0][0]})"
                    )
                elif self[P1][0] > n_cells:
                    raise ValueError(
                        f"I1 interface out of bounds: {self[P1][0]} > {n_cells}"
                    )
        return
