import os
import shutil

from pathlib import Path
from zipfile import ZipFile
from functools import wraps


def work_in_zipped_dir(zip_path: Path):
    """
    Extract some data from a compressed folder, change directories to it if
    required, run the function then, if required change directories back out
    and then delete the generated folder
    """

    if not str(zip_path).endswith('.zip'):
        raise ValueError(f"Cannot work in {zip_path}. Not a .zip directory")

    def func_decorator(func):

        @wraps(func)
        def wrapped_function(*args, **kwargs):

            cwd = os.getcwd()
            dir_path = Path(cwd, zip_path.stem)

            with ZipFile(zip_path, 'r') as zip_folder:
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
