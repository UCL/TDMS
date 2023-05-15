import os
import sys
from glob import glob
from pathlib import Path

from data.input_generation.generate_test_input import (
    generate_test_input,
    start_MatlabEngine_with_extra_paths,
)

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))

TESTS_TO_REGEN = sorted(glob(LOCATION_OF_THIS_FILE + "/config_files/*.yaml"))
N_TESTS_TO_REGEN = len(TESTS_TO_REGEN)


def regenerate_all() -> None:
    """Regenerate the input data for all tests corresponding to the config files present in this directory."""
    print(f"Found {N_TESTS_TO_REGEN} config files")
    # Use a single MATLAB session to avoid startup and shutdown
    engine = start_MatlabEngine_with_extra_paths()

    for i, test_case in enumerate(TESTS_TO_REGEN):
        # Fetch config file
        config_file_loc = Path(test_case)
        print(f"Regenerating ({i+1}/{N_TESTS_TO_REGEN}) {config_file_loc}")
        # Regenerate input data
        generate_test_input(config_file_loc, engine)

    # Quit the engine we started
    engine.quit()
    return


if __name__ == "__main__":
    regenerate_all()
    sys.exit(0)
