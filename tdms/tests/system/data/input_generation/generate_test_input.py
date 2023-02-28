import os
from glob import glob
from pathlib import Path
from typing import Literal

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
# Default values to use when an optional argument is not present in a config.yaml file
DEFAULT_VALUES = {
    "obstacle": 'fs',
    "obstacle_radius": 15.e-6,
}

def _create_temporary_filesetup(input_filename, temp_filesetup_name) -> None:
    """Generate the tempoaray file to be passed to iteratefdtd_matrix in filesetup mode when an illumination file has also been specified. The "filesetup" mode-file is essentially identical to it's counterpart, but needs the efname and hfname variables NOT to be present.
    As such, the optimal way to get around this is to have Python copy the input_file (which in this instance is what needs to be used to setup the illumination), and then remove the definition of the efname & hfname variables in the copy to create the "filesetup" input.
    """
    # Copy the input_file (illumination-input) line-by-line to a temporary location for the filesetup-input
    # Do not copy across the lines that define the efname and hfname variables
    with open(input_filename, "r") as illumination_input:
        with open(temp_filesetup_name, "w") as filesetup_input:
            for line in illumination_input:
                # Remove any potential whitespace padding from the beginning of the line
                stripped_line = line.lstrip()
                # Write line, provided efname or hfname are not defined on it
                if not (
                    ("efname=" in stripped_line) or ("hfname=" in stripped_line)
                ):
                    filesetup_input.write(line)
    # filesetup_input is now ready, and identical to illumination_input save in the definition of efname and hfname
    return

def run_bscan(
    test_directory: Path | str, input_filename: Path | str, engine: MatlabEngine, obstacle: Literal['fs', 'cyl', 'sph', 'sc'] = DEFAULT_VALUES["obstacle"], obstacle_radius: float = DEFAULT_VALUES["obstacle_radius"], illfile_required: bool = False) -> None:
    """Wrapper for running the run_bscan MATLAB function in the MATLAB engine provided.

    MatlabEngine cannot parse Path objects so file and directory paths must be cast to string when calling.

    The bscan/ and matlab/ directories are assumed to already be in the
    includepath of the engine instance, so that the run_bscan and supporting
    MATLAB files can be called.

    The obstacle radius is the circular face radius for cylinders (cyl), sphere radius for spheres (sph), and is ignored by freespace (fs) and point-source (sc) obstacles.

    :param input_filename: The path to the input file that defines the variables iteratefdtd_matrix reads in
    :param engine: The MatlabEngine instance to call run_bscan within. If not provided, a new session will be started and ended once the call is complete.
    :param obstacle: The obstacle that is present in the simulation.
    :param obstacle_radius: Radius of the spatial obstacle in microns.
    :param illfile_required: Flags whether run_bscan requires a call to iteratefdtd_matrix in illsetup mode as well as filesetup mode.
    """
    # In the event that illsetup is required for this run, generate the temporary name for the input file to be passed to iteratefdtd_matrix in filesetup mode
    # This only occurs when illumination files are required in input-data regeneration.
    illfile_extra_file = ""
    if illfile_required:
        # If this is _not_ an empty string, we need to generate the illumination file from the input file
        # It is essentially identical, but needs the efname and hfname variables to be:
        # present in the workspace/file when calling with 'illsetup'
        # absent from the workspace/file when calling with 'filesetup'
        # Otherwise, the "input" file to the illsetup and input file for the .mat creation are identical. As such, the optimal way to get around this is to have Python pass the input_file via 'illsetup', then copy this file, remove the definition of the efname & hfname variables, then run 'filesetp' mode. We can then cleanup our "extra" file that we created.
        illfile_extra_file = os.path.splitext(input_filename)[0] + "temp_filesetup_file____.m"
        _create_temporary_filesetup(input_filename, illfile_extra_file)
    # function [] = run_bscan(test_directory, input_filename, non_fs_obstacle, illfile_extra_file, obstacle_radius)
    engine.run_bscan(str(test_directory), str(input_filename), obstacle, illfile_extra_file, obstacle_radius, nargout=0)

    # Cleanup the illfile_extra_file, if it was created
    if illfile_required:
        os.remove(illfile_extra_file)
    return


def start_MatlabEngine_with_extra_paths(
    working_directory: str | Path | None = None,
) -> MatlabEngine:
    """Starts a new MatlabEngine and adds the bscan/ and matlab/ folders to its path, which are required to be in scope when regenerating the input data.

    :param working_directory: The working directory to start the MatlabEngine in. Should be an absolute path. Defaults to the working directory of the currently executing script if not passed.
    :returns: MatlabEngine instance with the additional bscan/ and matlab/ files on the MATLABPATH.
    """
    engine = matlab.start_matlab(MATLAB_STARTUP_OPTS)
    # Change to requested working directory if provided
    if working_directory:
        engine.cd(str(working_directory))
    # Append tdms scripts and functions to MATLABPATH
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
    # Fetch the spatial obstacles as a list
    obstacles = generation_info["spatial_obstacles"]
    # Fetch the non-freespace obstacle
    # Create a copy so we don't destroy the original list when finding the non-freespace obstacle
    non_freespace_obstacle = list(obstacles) 
    while "fs" in obstacles:
        obstacles.remove("fs")
    # Now obstacles should be left with only 1 element (the desired obstacle), if not then the non-freespace obstacle is non-unique and we error
    if len(obstacles) != 1:
        raise RuntimeError(
            f"Error: non-freespace obstacle is not unique ({obstacles})"
        )
    # If we didn't error, the only obstacle that remains is the non-freespace obstacle
    non_freespace_obstacle = non_freespace_obstacle[0]
    # Fetch the obstacle radius if it is present, otherwise use the default value
    if "obstacle_radius" in generation_info.keys() and (generation_info["obstacle_radius"] != None):
        obstacle_radius = float(generation_info["obstacle_radius"])
    else:
        obstacle_radius = DEFAULT_VALUES["obstacle_radius"]
    # Fetch whether illsetup mode is required
    if generation_info["illsetup"]:
        illsetup_requred = True
    else:
        # Cast things like None to bools, so typehints and behaviour is consistent
        illsetup_requred = False

    # Determine if we need to create our own MATLAB session
    # Explicit instance check since MatlabEngine may not have implicit casts/ interpretations
    engine_provided = isinstance(engine, MatlabEngine)
    if not engine_provided:
        # Start a new Matlab engine operating in the test directory
        engine = start_MatlabEngine_with_extra_paths(working_directory=test_dir)

    run_bscan(test_dir, input_file, engine, non_freespace_obstacle, obstacle_radius, illsetup_requred)

    # Quit our temporary MATLAB session, if we started one
    if not engine_provided:
        engine.quit()
    # Cleanup auxillary .mat files that are placed into this directory
    for aux_mat in sorted(glob(LOCATION_OF_THIS_FILE + "/*.mat")):
        os.remove(aux_mat)
        
    return