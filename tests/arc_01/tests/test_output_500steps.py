import os
import numpy as np
import h5py

from pathlib import Path
from utils import work_in_zipped_dir


THIS_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
OUT_DIR_PATH = Path(THIS_DIR_PATH, "..", "out")


@work_in_zipped_dir(Path(THIS_DIR_PATH, "test_data.zip"))
def test_output_as_expected():
    """
    For all the .mat files in this directory ensure that the generated file in
    the out/ directory matches. Where a match is defined in terms of closeness
    of all the arrays in the file
    """

    for filepath in OUT_DIR_PATH.iterdir():

        reference_filepath = Path(os.getcwd(), filepath.name)
        assert reference_filepath.exists()

        reference_file = HDF5File(reference_filepath)
        assert HDF5File(filepath).matches(reference_file)

    return None


class HDF5File(dict):
    """
    HDF5 file created from a .mat. Provides a dictionary of key value pairs
    where the keys are variable names defined in the .mat file and the values
    numpy arrays created from them.
    """

    COMPLEX_TYPE_STRING = r"[('real', '<f8'), ('imag', '<f8')]"

    def __init__(self, filepath: Path):
        super().__init__()

        with h5py.File(filepath, "r") as file:
            self.update({k: self.to_numpy_array(v) for k, v in file.items()})

    def to_numpy_array(self, dataset: h5py.Dataset) -> np.ndarray:
        """Convert a hdf5 dataset into a numpy array"""

        if str(dataset.dtype) == self.COMPLEX_TYPE_STRING:
            return self.tuple_to_complex(dataset)

        return np.array(dataset)

    @staticmethod
    def tuple_to_complex(dataset: h5py.Dataset) -> np.ndarray:
        """
        Convert an array who's last dimension is a tuple of reals into a complex
        valued array
        """

        array = np.array(dataset)
        shape = array.shape
        array = array.flatten()
        array = np.fromiter((complex(*xi) for xi in array), complex)

        return array.reshape(shape)

    def matches(self, other: "HDF5File", rtol=1e-4) -> bool:
        """
        Does this file match another. All arrays must be within a rtol
        to the other, where rtol is the relative difference between the two
        arrays
        """

        for key, value in self.items():

            if key not in other:
                return False  # Key did not match

            other_value = other[key]

            if value.shape != other_value.shape:
                return False  # Shapes did not match

            if relative_mean_squared_difference(value, other_value) > rtol:
                print(f"{key} was not within {rtol} to the reference (rel MSD)")
                return False

        return True


def relative_mean_squared_difference(a: np.ndarray, b: np.ndarray) -> float:
    """Calculate the relative mean square difference between two arrays"""

    mean_sq_a = np.mean(np.square(a))
    mean_sq_b = np.mean(np.square(b))

    if mean_sq_a < 1e-16 and mean_sq_b < 1e-16:
        return 0.0  # Prevent division by zero for tensors filled with zeros

    return np.mean(np.square(a - b)) / max(mean_sq_a, mean_sq_b)
