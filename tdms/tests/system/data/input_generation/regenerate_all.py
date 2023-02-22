import os
import sys
from glob import glob
from pathlib import Path

from generate_test_input import GenerationData

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))

TESTS_TO_REGEN = glob(LOCATION_OF_THIS_FILE + "/*.yaml")
N_TESTS_TO_REGEN = len(TESTS_TO_REGEN)


def regenerate_all() -> None:
    print(f"Found {N_TESTS_TO_REGEN} config files")
    for i, test_case in enumerate(TESTS_TO_REGEN):
        # Fetch config file
        config_file_loc = Path(test_case)
        print(f"Regenerating ({i+1}/{N_TESTS_TO_REGEN}) {config_file_loc}")
        # Regenerate input data
        regenerator = GenerationData(config_file_loc)
        regenerator.generate()
    return


if __name__ == "__main__":
    regenerate_all()
    sys.exit(0)
