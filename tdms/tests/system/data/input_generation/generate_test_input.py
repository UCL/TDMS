import os
from glob import glob
from pathlib import Path

import matlab.engine as matlab
import yaml
from matlab.engine import MatlabEngine

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))
# Additional options for running matlab on the command-line
MATLAB_OPTS_LIST = ["-nodisplay", "-nodesktop", "-nosplash", "-r"]
MATLAB_STARTUP_OPTS = " ".join(MATLAB_OPTS_LIST)
# Paths to matlab functions not in LOCATION_OF_THIS_FILE
MATLAB_EXTRA_PATHS = [
    os.path.abspath(LOCATION_OF_THIS_FILE + "/bscan"),
    os.path.abspath(LOCATION_OF_THIS_FILE + "/matlab"),
]


def run_bscan(
    test_directory: Path | str, input_filename: Path | str, engine: MatlabEngine
) -> None:
    """Wrapper for running the run_bscan MATLAB function in the MATLAB engine provided.

    MatlabEngine cannot parse Path objects so file and directory paths must be cast to string when calling.

    The bscan/ and matlab/ directories are assumed to already be in the
    includepath of the engine instance, so that the run_bscan and supporting
    MATLAB files can be called.
    """
    # function [] = run_bscan(test_directory, input_filename)
    engine.run_bscan(str(test_directory), str(input_filename), nargout=0)
    return


def start_MatlabEngine_with_extra_paths() -> MatlabEngine:
    """Starts a new MatlabEngine and adds the bscan/ and matlab/ folders to its path, which are required to be in scope when regenerating the input data."""
    engine = matlab.start_matlab(MATLAB_STARTUP_OPTS)
    for path in MATLAB_EXTRA_PATHS:
        engine.addpath(path)
    return engine


def generate_test_input(
    config_filepath: Path | str, engine: MatlabEngine | None = None
) -> None:
    """(re)Generates the input data (.mat files) contained in the config file, using the MATLAB session provided.

    This function is equivalent to running the run_{pstd,fdtd}_bscan.m scripts on the (test corresponding to the) config file in question.

    :param config_filepath: The path to the config file containing information about this system test
    :param engine: The MATLAB session to run the run_bscan function within. A session will be created and quit() if one is not provided.
    """
    with open(config_filepath, "r") as file:
        config_data = yaml.safe_load(file)

    # ID of the test we are generating input data for
    test_id = config_data["test_id"]
    # Absolute path to the directory into which the input data should be placed
    test_dir = Path(LOCATION_OF_THIS_FILE, "arc_" + test_id)
    # Ensure that the directory to place the output into exists, or create it otherwise
    if not test_dir.exists():
        print(f"The Path {test_dir} does not exist - creating now")
        os.mkdir(test_dir)
    elif not test_dir.is_dir():
        raise RuntimeError(f"{test_dir} is not a directory!")
    # else: the directory already exists, we don't need to do anything

    # Extract necessary input data generation information
    generation_info = config_data["input_generation"]
    # Fetch the location of the input file that generates the binary .mat input
    input_file = Path(LOCATION_OF_THIS_FILE, generation_info["input_file"])
    if not input_file.exists():
        raise RuntimeError(f"{input_file} does not exist")
    # Fetch the spatial obstacles
    obstacles = generation_info["spatial_obstacles"]

    # Determine if we need to create our own MATLAB session
    engine_provided = isinstance(engine, MatlabEngine)
    if not engine_provided:
        # Explicit instance check since MatlabEngine may not have implicit casts/ interpretations
        engine = start_MatlabEngine_with_extra_paths()

    run_bscan(test_dir, input_file, engine)

    # Quit our temporary MATLAB session, if we started one
    if not engine_provided:
        engine.quit()
    # cleanup auxillary .mat files that are placed into this directory
    for aux_mat in sorted(glob(LOCATION_OF_THIS_FILE + "/*.mat")):
        os.remove(aux_mat)
