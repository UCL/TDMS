"""
Common utilities for running TDMS system tests
"""
import os
import h5py
import shutil
import numpy as np

from dataclasses import dataclass
from urllib import request
from platform import system
from typing import Union
from pathlib import Path
from zipfile import ZipFile
from functools import wraps
from subprocess import Popen, PIPE

executable_name = "tdms.exe" if system() == "Windows" else "tdms"
executable_path = shutil.which(executable_name)

if Path(executable_name).is_file():
    # If the executable exists in the current working directory use that
    executable_path = str(Path(executable_name).absolute())


class HDF5File(dict):
    """
    HDF5 file created from a .mat. Provides a dictionary of key value pairs
    where the keys are variable names defined in the .mat file and the values
    numpy arrays created from them.
    """

    COMPLEX_TYPE_STRING = r"[('real', '<f8'), ('imag', '<f8')]"

    def __init__(self, filepath: Union[str, Path]):
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

    def matches(self, other: "HDF5File", rtol=1e-10) -> bool:
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
                print(f"Array shape in {key} was not the same. "
                      f"{value.shape} ≠ {other_value.shape}")
                return False  # Shapes did not match

            r_ms_diff = relative_mean_squared_difference(value, other_value)
            if r_ms_diff > rtol:
                print(f"{key} was not within {rtol} to the reference. "
                      f"relative MSD = {r_ms_diff:.8f})")
                return False

        return True


def relative_mean_squared_difference(a: np.ndarray, b: np.ndarray) -> float:
    """Calculate the relative mean square difference between two arrays"""

    mean_sq_a = np.mean(np.square(np.abs(a)))
    mean_sq_b = np.mean(np.square(np.abs(b)))

    if mean_sq_a < 1e-16 and mean_sq_b < 1e-16:
        return 0.0  # Prevent division by zero for tensors filled with zeros

    return np.mean(np.square(np.abs(a - b))) / max(mean_sq_a, mean_sq_b)


def work_in_zipped_dir(zip_path: Path):
    """
    Function decorator to extract some data from a compressed folder into a
    directory, run the function, change directories to the initial working
    directory and then delete the generated folder.
    """

    if not str(zip_path).endswith(".zip"):
        raise ValueError(f"Cannot work in {zip_path}. Not a .zip directory")

    def func_decorator(func):
        @wraps(func)
        def wrapped_function(*args, **kwargs):

            cwd = os.getcwd()
            dir_path = Path(cwd, zip_path.stem)

            with ZipFile(zip_path, "r") as zip_folder:
                zip_folder.extractall(dir_path)

            os.chdir(dir_path)

            try:
                result = func(*args, **kwargs)

            finally:
                os.chdir(cwd)
                shutil.rmtree(dir_path)

            return result

        return wrapped_function

    return func_decorator


@dataclass
class Result:

    return_code: int
    stdout: str


def run_tdms(*args) -> Result:
    """
    Run the tdms executable. Requires a tdms executable within the working
    directory or $PATH.
    """

    if executable_path is None:
        raise AssertionError("Failed to run tdms. Not found in either current "
                             "working directory or $PATH")

    p = Popen([executable_path, *args], stdout=PIPE)
    stdout, _ = p.communicate()

    return Result(p.returncode, stdout.decode())


def download_data(url: str, to: Path) -> None:
    """Download binary data present at a URL to a file"""
    print("Downloading: ", url)

    cwd = Path().cwd()
    if not to.parent.exists():
        to.parent.mkdir()

    os.chdir(to.parent)
    response = request.urlopen(url)

    with open(to, "wb") as file:
        file.write(response.read())

    os.chdir(cwd)

    return None
