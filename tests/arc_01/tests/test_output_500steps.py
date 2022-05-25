import os
import numpy as np
import h5py

from pathlib import Path


THIS_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
COMPLEX_TYPE_STRING = r"[('real', '<f8'), ('imag', '<f8')]"


def test_output_as_expected():
    raise NotImplementedError


class HDF5File:
    """
    HDF5 file built from a .mat. Provides a dictionary of key value pairs
    where the keys are variable names defined in the .mat file and the values
    numpy arrays created from them
    """

    def __init__(self, filepath: Path):

        with h5py.File(filepath, 'r') as file:
            self._dict = {k: self.to_numpy_array(v) for k, v in file.items()}

    def __getitem__(self, item):
        return self._dict[item]

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

    def to_numpy_array(self, dataset: h5py.Dataset) -> np.ndarray:
        """Convert a hdf5 dataset into a numpy array"""

        if str(dataset.dtype) == COMPLEX_TYPE_STRING:
            return self.tuple_to_complex(dataset)

        return np.array(dataset)


if __name__ == '__main__':

    path = Path(THIS_DIR_PATH, '..', 'out', 'pstd_fs.mat')
    f = HDF5File(path)
