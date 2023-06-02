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

    # Create a group under root
    read_in_test = file.require_group("read_in_test")

    # Create test data to read in
    consecutive_numbers = np.arange(0, 12, dtype=float)

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

    file.attrs["file_attribute"] = 1

    file.close()
    return


if __name__ == "__main__":
    create_hdf5_test_file()
    sys.exit(0)
