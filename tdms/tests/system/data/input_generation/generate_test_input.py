import os
from glob import glob
from io import StringIO
from pathlib import Path
from typing import Any, Tuple

import matlab.engine as matlab
import yaml
from matlab.engine import MatlabEngine

from .matfile_option_edit import edit_mat_file

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))
INPUT_M_FILE_LOCATION = LOCATION_OF_THIS_FILE + "/input_files/"
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
    "obstacle": "fs",
    "obstacle_radius": 15.0e-6,
    "illsetup": False,
    "refind": 1.42,
    "calc_tdfield": False,
}
OPTIONAL_ARGS = DEFAULT_VALUES.keys()


def bscan_options(output_name: str, generation_info: dict[str, Any]) -> dict[str, Any]:
    """Creates a dictionary that can be passed as a structure to run_bscan.m, containing the options for producing the .mat file that run_bscan expects. Default values are used if the generation information does not contain the corresponding optional fields.

    :param output_name: The name of the .mat file that is to be produced.
    :param generation_info: Information from a config.yaml file containing the options to use when generating the output_name.mat data file.
    :returns: A dictionary that can be passed to run_bscan.m as the "options" parameter, to generate the output_name.mat file.
    """
    # Assume the default values unless told otherwise
    options = dict(DEFAULT_VALUES)
    # Populate the optional arguments that were passed, ignoring empty fields
    for arg, passed_value in generation_info.items():
        if (arg in OPTIONAL_ARGS) and (passed_value != None):
            # We should cast to the correct type before storing the value, since MatlabEngine throws errors when it recieves the wrong datatypes
            expected_type = type(DEFAULT_VALUES[arg])
            options[arg] = expected_type(passed_value)
    # Populate the name of the .mat file that will be produced
    options["output_name"] = output_name

    return options


def run_bscan(
    test_directory: Path | str,
    matfile_to_produce: str,
    generation_info: dict[str, Any],
    engine: MatlabEngine,
) -> Tuple[list[str], StringIO, StringIO]:
    """Wrapper for running the run_bscan MATLAB function in the MATLAB engine provided.

    MatlabEngine cannot parse Path objects so file and directory paths must be cast to string when calling.

    The bscan/ and matlab/ directories are assumed to already be in the
    includepath of the engine instance, so that the run_bscan and supporting
    MATLAB files can be called.

    The obstacle radius is the circular face radius for cylinders (cyl), sphere radius for spheres (sph), and is ignored by freespace (fs) and point-source (sc) obstacles.

    The first output returned is the output of the run_bscan.m function, as a list. Additional outputs are the stdout and stderr of the MATLAB call.

    :param input_filename: The path to the input file that defines the variables iteratefdtd_matrix reads in
    :param engine: The MatlabEngine instance to call run_bscan within. If not provided, a new session will be started and ended once the call is complete.
    :param obstacle: The obstacle that is present in the simulation.
    :param obstacle_radius: Radius of the spatial obstacle in microns.
    :param illsetup: Flags whether run_bscan requires a call to iteratefdtd_matrix in illsetup mode as well as filesetup mode.
    :param calc_tdfield: Flags whether run_bscan must setup a time-domain field and pass the resulting .mat file into iteratefdtd_matrix.
    :returns: (stdout, stderr) produced by the MatlabEngine whilst running.
    """
    # The .m input file that will be read by iteratefdtd_matrix
    input_filename = os.path.abspath(
        INPUT_M_FILE_LOCATION + generation_info["input_file"]
    )
    # Create the options dictionary (which will be converted to a struct) to pass to run_bscan.m
    options = bscan_options(matfile_to_produce, generation_info)
    # Create IO objects to capture stdout and stderr from MatlabEngine, to avoid polluting the terminal
    matlab_stdout = StringIO()
    matlab_stderr = StringIO()

    # function [freespace_output_file, obstacle_output_file] = run_bscan(test_directory, input_filename, non_fs_obstacle, illfile_extra_file, obstacle_radius, calc_tdfield)
    engine.run_bscan(
        str(test_directory),
        str(input_filename),
        options,
        nargout=0,
        stdout=matlab_stdout,
        stderr=matlab_stderr,
    )

    # Return stdout and stderr messages
    return matlab_stdout, matlab_stderr


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

    # Names of the .mat files to produce
    mats_to_produce = [key for key in config_data.keys() if key != "test_id"]
    # Some or all entries in mats_to_produce might be missing the .mat extension. Keeping the original key names from the config.yaml file is necessary for input generation purposes (else keyerrors are thrown) and MATLAB will auto-append the extension on saving, so there will be no name mismatches.
    # However we need to make sure all our output files have the .mat extension provided, and are absolute paths, when we clean up after producing the outputs, otherwise they will be deleted along with other .mat artifacts.
    mats_produced_with_extension = [str(test_dir / mfile) for mfile in mats_to_produce]
    for i, file in enumerate(mats_produced_with_extension):
        if file[-4:] != ".mat":
            mats_produced_with_extension[i] += ".mat"

    # Determine if we need to create our own MATLAB session
    # Explicit instance check since MatlabEngine may not have implicit casts/ interpretations
    engine_provided = isinstance(engine, MatlabEngine)
    if not engine_provided:
        # Start a new Matlab engine operating in the test directory
        engine = start_MatlabEngine_with_extra_paths(working_directory=test_dir)
    # Capture the working directory of the engine (for cleanup)
    matlab_working_directory = engine.pwd(nargout=1)

    # Loop over all input .mat files to be produced by this config file
    # Input files that lack the "adjust" key require calls to run_bscan, and must be produced first.
    bscan_matfiles = list(mats_to_produce)
    adjust_matfiles = list(mats_to_produce)
    for matfile in mats_to_produce:
        if (
            "adjust" in config_data[matfile].keys()
            and config_data[matfile]["adjust"] != None
        ):
            bscan_matfiles.remove(matfile)
        else:
            adjust_matfiles.remove(matfile)
    # Those that possess the "adjust" key need to be done afterwards, to ensure the file they depend on is present.
    for matfile in bscan_matfiles:
        _, _ = run_bscan(test_dir, matfile, config_data[matfile], engine)
    for matfile in adjust_matfiles:
        edit_mat_file(test_dir, matfile, config_data[matfile])

    # Quit our temporary MATLAB session, if we started one
    if not engine_provided:
        engine.quit()
    # Cleanup auxillary .mat files that are placed into this directory, and the MATLAB working directory
    auxillary_matfiles = (
        set(glob(str(test_dir) + "/*.mat"))
        | set(glob(matlab_working_directory + "/*.mat"))
    ) - set(mats_produced_with_extension)
    for aux_mat in auxillary_matfiles:
        os.remove(aux_mat)

    return
