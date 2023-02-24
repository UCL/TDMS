import os
import sys
from glob import glob
from pathlib import Path
from typing import Union

from generate_test_input import GenerationData

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))

TESTS_TO_REGEN = glob(LOCATION_OF_THIS_FILE + "/*.yaml")
N_TESTS_TO_REGEN = len(TESTS_TO_REGEN)


def regenerate_test(config_file_location: Union[Path, str]) -> tuple[int, str]:
    """Given the config file of a particular test, regenerate the input data for that test."""
    regenerator = GenerationData(config_file_location)
    return_code = regenerator.generate()
    return return_code, regenerator.matlab_command.stdout


def regenerate_all() -> None:
    """Regenerate the input data for all tests corresponding to the config files present in this directory."""
    print(f"Found {N_TESTS_TO_REGEN} config files")
    for i, test_case in enumerate(TESTS_TO_REGEN):
        # Fetch config file
        config_file_loc = Path(test_case)
        print(f"Regenerating ({i+1}/{N_TESTS_TO_REGEN}) {config_file_loc}")
        # Regenerate input data
        regenerate_test(config_file_loc)
    return


if __name__ == "__main__":
    regenerate_all()
    sys.exit(0)