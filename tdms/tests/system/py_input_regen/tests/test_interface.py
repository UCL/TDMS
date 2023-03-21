import numpy as np
import pytest

from ..interface import Interface


def test_keys() -> None:
    """Check that:
    Valid plane keys can be passed in either upper or lowercase.
    Invalid keys throw an error.
    """
    i = Interface(i0=[1, 1], k1=[5, 1])
    # The keys we provided should be in the object despite being passed as lowercase.
    assert {"I0", "K0"}.issubset(i.keys())
    # If we try to create / assign a non-plane item we should recieve an error
    with pytest.raises(KeyError):
        i["foo"] = [3, 0]

    return


def test_items() -> None:
    """Check that:
    Items can only be created if they can be cast to (2,) int-numpy arrays.
    Missing planes are automatically added to the object upon construction, with zero boundary conditions.
    """
    # Create an interface and leave out some keys
    i = Interface(K0=[1, 1], J1=[3, 0])
    # The object should still have 6 keys even though we only assigned 2 of them
    assert len((i.keys())) == 6
    # Assert that the additional planes are all zero-boundary condition interfaces
    for plane in ["i0", "i1", "j0", "k1"]:
        assert not i.non_zero_interface(
            plane
        ), f"{plane} is a non-zero interface, but should've been initialised as such."
    # Assert that we cannot assign items that cannot be cast to (2,) int arrays
    with pytest.raises(ValueError):
        i["K0"] = [1, 2, 3]
    with pytest.raises(ValueError):
        i["K0"] = [3]
    with pytest.raises(ValueError):
        i["K0"] = ["flibble", 1]

    return


def test_validation_and_cell_ranges() -> None:
    """Check that:
    Non-zero planes must satisfy the inequality 1 <= P0[0] <= P1[0] <= P.
    When valid, we can extract the cell_range and change_coords appropriately
    """
    i = Interface(k0=[1, 1], k1=[10, 1], j0=[0, 0], j1=[4, 1])
    # This validation will pass, since both J0 and J1 are not non-zero planes and so are ignored.
    i.validate_against_n_cells(K=10, J=2)
    # We should be able to pull out the cell range over K too
    assert np.all(i.cell_range("K") == np.arange(1, 10 + 1, step=1, dtype=int))
    # However, we shouldn't be able to extract the cell range over J
    with pytest.raises(ValueError):
        i.cell_range("J")

    # Now, "activate" the J plane and check for inequality errors
    with pytest.raises(ValueError):
        i["j0"] = [0, 1]
        i.validate_against_n_cells(J=5)
    with pytest.raises(ValueError):
        i["j0"] = [2, 1]
        i["j1"] = [1, 1]
        i.validate_against_n_cells(J=5)
    with pytest.raises(ValueError):
        i["j0"] = [1, 1]
        i["j1"] = [4, 1]
        i.validate_against_n_cells(J=3)

    return


def test_change_coords() -> None:
    """Check that we can change indexing schemes correctly."""
    # Create "fake" PML
    PML = {
        "Dxl": 1,
        "Dxu": 2,
        "Dyl": 3,
        "Dyu": 4,
        "Dzl": 5,
        "Dzu": 6,
    }

    # Create an interface in local coordinates
    originally_local = Interface(K0=[1, 1])
    assert originally_local.in_coords == "local"

    # Change coords, so now we are in global coords
    originally_local.change_coords("global", PML)
    assert originally_local.in_coords == "global"
    assert originally_local["K0"][0] == 6  # 1+Dzl=6
    assert originally_local["J0"][0] == PML["Dyl"]

    # Transform back to local coordinates
    originally_local.change_coords("local", PML)
    assert originally_local.in_coords == "local"
    assert originally_local["K0"][0] == 1
    assert originally_local["J0"][0] == 0

    return
