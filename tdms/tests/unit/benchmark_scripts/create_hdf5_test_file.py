import os
import sys

import h5py
import numpy as np

FNAME_TO_CREATE = os.path.abspath(
    os.path.dirname(__file__) + "/unit_test_data/hdf5_test_file.hdf5"
)


def create_hdf5_test_file() -> None:
    """ """
    # Create the directory for the file (if it doesn't exist), and the file itself
    if not os.path.exists(os.path.dirname(FNAME_TO_CREATE)):
        os.mkdir(os.path.dirname(FNAME_TO_CREATE))
    file = h5py.File(FNAME_TO_CREATE, "w")

    # Create test data to read in
    consecutive_numbers = np.arange(0, 12, dtype=float)

    # Create a group under root
    read_in_test = file.require_group("read_in_test")

    # Populate group with test data
    read_in_test.create_dataset(
        "vector_int", data=consecutive_numbers, shape=(12,), dtype=int
    )
    read_in_test.create_dataset(
        "matrix_double", data=consecutive_numbers, shape=(2, 6), dtype=float
    )
    read_in_test.create_dataset(
        "tensor_double", data=consecutive_numbers, shape=(2, 3, 2), dtype=float
    )

    # Create data for an XYZVector
    x_array = np.array([0.1, 0.2, 0.3])
    y_array = np.array([0.4, 0.5, 0.6])
    z_array = np.array([0.7, 0.8, 0.9])

    xyz_vector_group = read_in_test.create_group("XYZVector")
    xyz_vector_group.create_dataset("xyz_x", data=x_array)
    xyz_vector_group.create_dataset("xyz_y", data=y_array)
    xyz_vector_group.create_dataset("xyz_z", data=z_array)

    # Create & populate the group that mimics MATLAB empty arrays
    # Deliberately include some data here to stress that emptiness is based off the presence of an attribute
    flag_as_empty = file.require_group("flag_as_empty")
    flag_as_empty.attrs["MATLAB_empty"] = np.array(1, dtype=np.uint8)
    flag_as_empty.create_dataset(
        "consecutive_numbers", data=consecutive_numbers[0:3], shape=(3,), dtype=float
    )

    # Create a top-level dataset that will be empty, and give it the MATLAB_empty attribute.
    # But set this to be false, to emphasise dependence on this value.
    not_marked_empty = file.create_dataset(
        "not_marked_empty", data=consecutive_numbers[0:0], shape=(0,), dtype=float
    )
    not_marked_empty.attrs["MATLAB_empty"] = np.array(0, dtype=np.uint8)

    # Create an attribute at the root of the file for debugging and testing purposes
    file.attrs["file_attribute"] = 1

    file.close()
    return


if __name__ == "__main__":
    create_hdf5_test_file()
    sys.exit(0)
