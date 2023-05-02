import os
from glob import glob
from io import StringIO
from pathlib import Path
from typing import Any, Tuple

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
    "obstacle": "fs",
    "obstacle_radius": 15.0e-6,
    "illsetup": False,
    "ill_filesetup": ".__TDMS_input_regen_extra_filesetup_file.m",
    "refind": 1.42,
    "calc_tdfield": False,
}
OPTIONAL_ARGS = DEFAULT_VALUES.keys()


def _create_temporary_filesetup(file_to_remove_eh_name, file_to_place_into) -> None:
    """Generates a temporary file that will be used as input to iteratefdtd_matrix in filesetup mode, along with an illumination file.

    The file in filesetup mode is almost identical to the original input file, but with empty strings set for the efname and hfname variables. To achieve this, the input file (used to set up the illumination) is copied by Python, and the efname and hfname variables are modified to create the file in filesetup mode.

    :param file_to_remove_eh_name: The name of the file to copy the lines of, removing efname and hfname.
    :param file_to_place_into: The name of the file to write to.
    """
    # Copy the input_file (illumination-input) line-by-line to a temporary location for the filesetup-input
    # Do not copy across the lines that define the efname and hfname variables
    with open(file_to_remove_eh_name, "r") as illumination_input:
        with open(file_to_place_into, "w") as filesetup_input:
            for line in illumination_input:
                # Remove any potential whitespace padding from the line
                # This avoids funny business if there's whitespace around the = symbol where {ef,hf}name are defined
                stripped_line = line.replace(" ", "")
                # Write line, provided efname or hfname are not defined on it
                if ("efname=" not in stripped_line) and (
                    "hfname=" not in stripped_line
                ):
                    filesetup_input.write(line)
                elif "efname=" in stripped_line:
                    filesetup_input.write("efname = '';\n")
                elif "hfname=" in stripped_line:
                    filesetup_input.write("hfname = '';\n")
    # filesetup_input is now ready, and identical to illumination_input save in the definition of efname and hfname
    return


def bscan_options(output_name: str, generation_info: dict[str, Any]) -> dict[str, Any]:
    """Creates a dictionary that can be passed as a structure to run_bscan.m, containing the options for producing the .mat file that run_bscan expects. Default values are used if the generation information does not contain the corresponding optional fields.

    :param output_name: The name of the .mat file that is to be produced.
    :param generation_info: Information from a config.yaml file containing the options to use when generating the output_name.mat data file.
    :returns: A dictionary that can be passed to run_bscan.m as the "options" parameter, to generate the output_name.mat file.
    """
    # Assume the default values unless told otherwise
    options = dict(DEFAULT_VALUES)
    # Populate the optional arguments that were passed
    for arg, passed_value in generation_info.items():
        if arg in OPTIONAL_ARGS:
            options[arg] = passed_value
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
    input_filename = generation_info["input_file"]
    # Create the options dictionary (which will be converted to a struct) to pass to run_bscan.m
    options = bscan_options(matfile_to_produce, generation_info)
    # In the event that illsetup is required for this run, generate the temporary name for the input file to be passed to iteratefdtd_matrix in filesetup mode
    # This only occurs when illumination files are required in input-data regeneration.
    if options["illsetup"]:
        # We need to generate the illumination file from the input file.
        # It is essentially identical, but needs the efname and hfname variables to be:
        # present in the workspace/file when calling with 'illsetup'
        # absent from the workspace/file when calling with 'filesetup'
        # Otherwise, the "input" file to the illsetup and input file for the .mat creation are identical. As such, the optimal way to get around this is to have Python pass the input_file via 'illsetup', then copy this file, remove the definition of the efname & hfname variables, then run 'filesetp' mode. We can then cleanup our "extra" file that we created.
        _create_temporary_filesetup(input_filename, options["ill_filesetup"])
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

    # Cleanup the illfile_extra_file, if it was created
    if options["illsetup"]:
        os.remove(options["ill_filesetup"])
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
    # Names of the .mat files to produce
    mats_to_produce = [key for key in config_data.keys() if key != "test_id"]
    # Append .mat extension to the filenames of the inputs to be created, if they are not there already
    for i, file in enumerate(mats_to_produce):
        if file[-4:] != ".mat":
            mats_to_produce[i] += ".mat"
    # Absolute path to the directory into which the input data should be placed
    test_dir = Path(LOCATION_OF_THIS_FILE, "arc_" + test_id)
    # Ensure that the directory to place the output into exists, or create it otherwise
    if not test_dir.exists():
        print(f"The Path {test_dir} does not exist - creating now")
        os.mkdir(test_dir)
    elif not test_dir.is_dir():
        raise RuntimeError(f"{test_dir} is not a directory!")
    # else: the directory already exists, we don't need to do anything

    # Determine if we need to create our own MATLAB session
    # Explicit instance check since MatlabEngine may not have implicit casts/ interpretations
    engine_provided = isinstance(engine, MatlabEngine)
    if not engine_provided:
        # Start a new Matlab engine operating in the test directory
        engine = start_MatlabEngine_with_extra_paths(working_directory=test_dir)
    # Capture the working directory of the engine (for cleanup)
    matlab_working_directory = engine.pwd(nargout=1)

    # Loop over all input .mat files to be produced by this config file
    for matfile in mats_to_produce:
        _, _ = run_bscan(test_dir, matfile, config_data[matfile], engine)

    # Quit our temporary MATLAB session, if we started one
    if not engine_provided:
        engine.quit()
    # Cleanup auxillary .mat files that are placed into this directory, and the MATLAB working directory
    auxillary_matfiles = (
        set(glob(test_dir + "/*.mat")) | set(glob(matlab_working_directory + "/*.mat"))
    ) - set(mats_to_produce)
    for aux_mat in auxillary_matfiles:
        os.remove(aux_mat)

    return
