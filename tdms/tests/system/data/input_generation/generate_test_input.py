import os
from glob import glob
from pathlib import Path
from typing import Union

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


engine: MatlabEngine = matlab.start_matlab(MATLAB_STARTUP_OPTS)


class BScanArguments:
    """Class to tidily contain the input arguments that need to be passed into the MATLAB function of the same name, to generate the input data.

    The bscan/run_bscan.m file contains the matlab function which generates the input data. Regrettably, we need to specify particular inputs to this script for each test, which requires us to translate the argument values as read from the config.yaml file into a long string of values in the correct order, which can in turn be called from matlab.

    BScanArguments is essentially a glorified dictionary, it's members sharing the names of the input arguments to run_bscan. It's create_bscan_argument can be used to convert the values that need to be passed into a string of the form:
    run_bscan(input_arguments_in_the_correct_order).
    """

    # The directory of the test whose input data is being generated
    test_directory: str
    # The path to the input file that defines the variables iteratefdtd_matrix.m takes in
    input_filename: str
    # The obstacles present in the simulations we are generating inputs for
    obstacles: list[str]
    # If the BScan call requires an illumination file to setup, provide a non-empty string here with the filename
    # Otherwise, assume that no illumination setup is necessary
    illfile_req: bool

    def __init__(
        self,
        test_directory: Union[Path, str],
        input_filename: Union[Path, str],
        obstacles: list[str],
        illfile: bool = False,
    ) -> None:
        """Initialise by unpacking values, and assigning defaults if necessary."""
        # Cast possible Path() to str() to guard against unsupported datatypes
        self.test_directory = str(test_directory)
        self.input_filename = str(input_filename)
        self.obstacles = obstacles
        self.illfile_req = illfile

        return

    def _temp_filesetup_name(self) -> str:
        """In the event that illsetup is required for this run, generate the temporary name for the input file to be passed to iteratefdtd_matrix in filesetup mode.

        Only occurs when illumination files are required in input-data regeneration. In this case, self.input_filename is passed to iteratefdtd_matrix in 'illsetup' mode.
        """
        # Input file name without the extension (.m)
        raw_input_name = os.path.splitext(self.input_filename)[0]
        # Append to filename an suitably unique extension
        raw_input_name += "temp_filesetup_file____.m"
        return raw_input_name

    def _create_temp_filesetup(self) -> None:
        """Generate the tempoaray file to be passed to iteratefdtd_matrix in filesetup mode when an illumination file has also been specified. The "filesetup" mode-file is essentially identical to it's counterpart, but needs the efname and hfname variables NOT to be present.

        As such, the optimal way to get around this is to have Python copy the input_file (which in this instance is what needs to be used to setup the illumination), and then remove the definition of the efname & hfname variables in the copy to create the "filesetup" input.
        """
        # Copy the input_file (illumination-input) line-by-line to a temporary location for the filesetup-input
        # Do not copy across the lines that define the efname and hfname variables
        with open(self.input_filename, "r") as illumination_input:
            with open(self._temp_filesetup_name(), "w") as filesetup_input:
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

    def non_freespace_obstacle(self) -> str:
        """Returns the shape of the scattering object present in the non-freespace simulation, as a string.

        Throws an error in the event that there is no such object, or if there are multiple such objects.
        """
        # Error if self.obstacles is not set
        if not self.obstacles:
            raise RuntimeError(f"Error: no obstacles specified in call to BScan!")
        # Create a copy of obstacles so we don't overwrite our member variable
        obstacles = list(self.obstacles)
        # We should be able to pop out 'fs' now - if we can't then we should error
        # We can use while here in case there are duplications of 'fs', as we self-update each time
        while "fs" in obstacles:
            obstacles.remove("fs")
        # Now obstacles should be left with only 1 element (the desired obstacle), if not then the non-freespace obstacle is non-unique and we error
        if len(obstacles) != 1:
            raise RuntimeError(
                f"Error: non-freespace obstacle is not unique ({obstacles})"
            )
        # If we didn't error, we can return the only element of obstacles, which is the non-fs obstacle
        return obstacles[0]

    def run_bscan(self, engine: MatlabEngine) -> None:
        """Runs the run_bscan function in the MatlabEngine provided.

        The bscan/ and matlab/ directories are assumed to already be in the includepath of the engine instance, so that the run_bscan and supporting matlab files can be called.
        """
        # If illumination file is needed, we need to copy the input file and remove the lines that define the efname and hfname as necessary
        illfile_extra_file = ""
        if self.illfile_req:
            # If this is _not_ an empty string, we need to generate the illumination file from the input file
            # It is essentially identical, but needs the efname and hfname variables to be:
            # present in the workspace/file when calling with 'illsetup'
            # absent from the workspace/file when calling with 'filesetup'
            # Otherwise, the "input" file to the illsetup and input file for the .mat creation are identical. As such, the optimal way to get around this is to have Python pass the input_file via 'illsetup', then copy this file, remove the definition of the efname & hfname variables, then run 'filesetp' mode. We can then cleanup our "extra" file that we created.
            illfile_extra_file = self._temp_filesetup_name()
            self._create_temp_filesetup()

        # function [] = run_bscan(test_directory, input_filename, non_fs_obstacle, illfile_extra_file)
        # pass nargout=0 to indicate no value is to be returned to Python
        engine.run_bscan(
            self.test_directory,
            self.input_filename,
            self.non_freespace_obstacle(),
            illfile_extra_file,
            nargout=0,
        )

        # Cleanup the temporary illumination file, if we created it
        if self.illfile_req:
            os.remove(illfile_extra_file)
        return


class MATLABEngineWrapper:
    """When we regenerate input data, we always need to add the bscan/ and matlab/ directories to the MATLAB instance's search path. correspondingly, we always want to kill the MATLAB instance after running the run_bscan function and generating the data.

    This class serves as a wrapper for that purpose. It stores instance(s) of the BScanArguments class, which it will run in sequence between the aforementioned addpath() setup and then engine shutdown. The .run() method performs exactly this.
    """

    # List containing the subsequent matlab commands to be executed by the interpreter
    bscan_calls: list[BScanArguments]

    # The MATLAB session that will run
    engine: MatlabEngine

    def __init__(self, bscans: Union[BScanArguments, list[BScanArguments]]) -> None:
        """Initialise with the provided (list of) BScanArguments."""
        self.bscan_calls = []
        # add the additional commands here
        if isinstance(bscans, BScanArguments):
            # handed BScanArgument, construct command from this
            self.bscan_calls = [bscans]
        elif isinstance(bscans, list):
            # validate this is a list of BScanArguments and not some other datatype
            non_bscan_items = sum([isinstance(b, BScanArguments) for b in bscans])
            if non_bscan_items > 0:
                raise RuntimeError(
                    f"Error: not all inputs are BScan calls ({non_bscan_items}/{len(bscans)})"
                )
            else:
                self.bscan_calls = bscans
        else:
            # I don't know what we've been handed, but it's not what I was expecting
            raise RuntimeError(
                f"Error: expected BScanArguments or list[BScanArguments] but got {type(bscans)}"
            )
        # do not start the engine yet
        self.engine = None
        return

    def _addpath_commands(self) -> None:
        """Adds the bscan/ and matlab/ directories to the session's path."""
        for path in MATLAB_EXTRA_PATHS:
            self.engine.addpath(path)
        return

    def _stop_engine(self) -> None:
        """Kills engine if it is still running."""
        self.engine.quit()
        return

    def run(self, kill_on_complete: bool = True) -> None:
        """Run the bscan arguments saved to the instance, in the same session.

        Engine is terminated if kill_on_complete is True, otherwise the engine is left running and manual cleanup is needed. It can be useful for debugging to leave the engine running, however it is recommended to stop it after we have finished using it.
        """
        # If we have no bscan arguments to run, don't bother starting the engine and report this
        if len(self.bscan_calls) == 0:
            raise RuntimeWarning(
                "No bscan calls specified in this instance. Engine not started."
            )
        else:
            # Start the engine
            self.engine = matlab.start_matlab(MATLAB_STARTUP_OPTS)
            # Add necessary paths
            self._addpath_commands()

            # Run every bscan call
            for b in self.bscan_calls:
                b.run_bscan(self.engine)

            # If requested, kill the instance
            if kill_on_complete:
                self._stop_engine()
        return


class GenerationData:
    """Handles the information and processes related to creating the input data for a particular test case.

    Information is read from a config.yaml file, and the contents of the "generation" field are validated against what is expected.

    From this information, the call to run_bscan is created. This command is the equivalent of running the old run_{pstd,fdtd}_bscan.m scripts on the test in question, however it is now sufficiently generalised so that we do not need an individual run file for each test.

    The input data for the particular test can then be generated using the .generate() method on the class instance.
    """

    # Location of the config file in case we need to go back to it
    _config_file_location: Path
    # Full option list from the config yaml file
    _generation_options: dict
    # The required options that the yaml file's generation field needs to contain
    _required_config_options: list[str] = ["input_file", "spatial_obstacles"]

    # ID of the test we are generating input data for
    test_id: str
    # Absolute path to the directory into which the input data should be placed
    test_dir: Path

    # Path to the input file that input data generation requires
    input_file: Path
    # The spatial obsticles that are present in the runs of this test
    obstacles: list[str]
    # Whether or not we need to run iteratefdtd_matrix in illsetup mode as well as filesetup mode
    illsetup_required: bool

    # The matlab instance that generates the input data for this test
    matlab_instance: MATLABEngineWrapper

    def __init__(self, config_file_location: Union[str, Path]) -> None:
        """Initialise instance by reading information from a config file"""
        # Fetch config file and options
        self._config_file_location = Path(config_file_location)
        with open(config_file_location, "r") as config_file:
            config = yaml.safe_load(config_file)

        # Get the ID of the test & thus the directory to save to
        self.test_id = config["test_id"]
        self.test_dir = Path(LOCATION_OF_THIS_FILE, ("arc_" + self.test_id))
        # Information about how this test's input data is generated
        self._generation_options = config["input_generation"]
        self._validate_config_options()

        # Setup the object based on the inputs
        self._setup_member_vars_from_config()

        # Setup the command from the member variables
        self.matlab_instance = self._setup_matlab_instance()
        return

    def _validate_config_options(self) -> None:
        """Validate that the config file contains the required fields for creation of the input data"""
        # Check for compulsory fields
        for required_option in self._required_config_options:
            if required_option not in self._generation_options.keys():
                raise RuntimeError(
                    f"{required_option} is not present in config file {self._config_file_location}"
                )
        return

    def _setup_member_vars_from_config(self) -> None:
        """Setup local variables by reading from the config file.

        These variables typically serve a purpose outside the generation process itself (EG specifying directories, locations to write to, etc).
        """
        # Fetch the location of the input file that generates the binary .mat input
        self.input_file = Path(
            LOCATION_OF_THIS_FILE, self._generation_options["input_file"]
        )
        if not self.input_file.exists():
            raise RuntimeError(f"{self.input_file} does not exist")
        # fetch the spatial obstacles
        self.obstacles = self._generation_options["spatial_obstacles"]
        # determine if illsetup is required in addition to filesetup
        if self._generation_options["illsetup"]:
            self.illsetup_required = True
        else:
            # Cast things like None to bools, so typehints and behaviour is consistent
            self.illsetup_required = False
        return

    def _find_or_create_test_dir(self) -> None:
        """Check if the directory to place the test data into already exists.
        If it does not, create it.
        """
        # check if the directory with the test data already exists
        if not self.test_dir.exists():
            print(f"The Path {self.test_dir} does not exist - creating now")
            os.mkdir(self.test_dir)
        elif not self.test_dir.is_dir():
            raise RuntimeError(f"{self.test_dir} is not a directory!")
        # else, the directory already exists, we don't need to do anything
        return

    def _setup_matlab_instance(self) -> MATLABEngineWrapper:
        """Setup the matlab command that will generate the input data"""
        bscan = BScanArguments(
            self.test_dir, self.input_file, self.obstacles, self.illsetup_required
        )
        return MATLABEngineWrapper(bscan)

    def generate(self) -> None:
        """Generate the input data to the test, as specified by this instance's member values.

        Returns the exit code of the sequence of matlab commands that were run."""
        # ensure that the directory to place the output into exists, or create it otherwise
        self._find_or_create_test_dir()

        # generate the input data for this test
        self.matlab_instance.run()

        # cleanup auxillary .mat files that are placed into this directory
        for aux_mat in sorted(glob(LOCATION_OF_THIS_FILE + "/*.mat")):
            os.remove(aux_mat)
        return
