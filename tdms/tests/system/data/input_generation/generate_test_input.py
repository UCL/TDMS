import os
from glob import glob
from pathlib import Path
from subprocess import PIPE, Popen, run
from typing import Union

import yaml

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))
# MATLAB executable path - need to locate this...
MATLAB_PATH = run(["which", "matlab"], stdout=PIPE).stdout.decode().removesuffix("\n")
# Additional options for running matlab on the command-line
MATLAB_OPTS = ["-nodisplay", "-nodesktop", "-nosplash", "-r"]
# Paths to matlab functions not in LOCATION_OF_THIS_FILE
MATLAB_EXTRA_PATHS = [
    os.path.abspath(LOCATION_OF_THIS_FILE + "/bscan"),
    os.path.abspath(LOCATION_OF_THIS_FILE + "/matlab"),
]


class BScanArguments:
    """Class to tidily contain the input arguments that need to be passed into the MATLAB function of the same name, to generate the input data.

    The bscan/run_bscan.m file contains the matlab function which generates the input data. Regrettably, we need to specify particular inputs to this script for each test, which requires us to translate the argument values as read from the config.yaml file into a long string of values in the correct order, which can in turn be called from matlab.

    BScanArguments is essentially a glorified dictionary, it's members sharing the names of the input arguments to run_bscan. It's create_bscan_argument can be used to convert the values that need to be passed into a string of the form:
    run_bscan(input_arguments_in_the_correct_order).
    """

    # The directory of the test whose input data is being generated
    test_directory: str
    # The path to the input file that defines the variables iteratefdtd_matrix.m takes in
    input_filename: Union[str, Path]

    def __init__(self, test_directory, input_filename) -> None:
        """Initialise by unpacking values."""
        self.test_directory = test_directory
        self.input_filename = input_filename
        return

    def create_bscan_argument(self) -> str:
        """Create the syntax of the call to run_bscan.m, but inserting the input-argument values in place.
        Also inserts "default" values if appropriate, to avoid doing this in matlab with varargin.

        The result is a string of the form:
        run_bscan(input_arguments_in_the_correct_order),

        where the input arguments are replaced with their values that need to be passed into the corresponding matlab call.
        """
        bscan_args = []
        bscan_args += ["'" + str(self.test_directory) + "'"]
        bscan_args += ["'" + str(self.input_filename) + "'"]

        bscan_arg_string = "run_bscan(" + ",".join(bscan_args) + ")"
        return bscan_arg_string


class MATLABCommand:
    """A class that allows us to create a matlab command that can be run in a shell.

    Serves as a wrapper for the input data generation. An instance of this class stores the commands we would like to execute in matlab_commands (list) to create the input data, in the desired order of execution. Commands for adding the paths to the .m files needed for data generation and to exit matlab after concluding are always the first and last elements of this list, and are automatically (pre/a)-ppended.

    The .run() method creates a shell instance that runs the MATLAB commands specified by the class instance, by converting the matlab_commands list into a string which is then passed to matlab on the command-line. The member variables stdout and return_code provide debugging information for the shell process that runs.
    """

    # List containing the subsequent matlab commands to be executed by the interpreter
    matlab_commands: list[str]
    # Stdout from running this command
    stdout: str
    # Shell exit code after running this matlab command
    return_code: int

    def __init__(self, matlab_args: Union[BScanArguments, list[str]]) -> None:
        """Initialise by being told the arguments that the run_bscan function will need.

        Then wrap a call to this command between addpath(), so that matlab can find the run_bscan function itself, and exit, so that matlab does not hang after running the bscan function.

        If the input is a BScanArguments object, the constructor can assemble the matlab call to run_bscan by itself.
        If the input is a list[str], then the constructor interprets this as a list of matlab commands to be run, in the list order.
        """
        self.matlab_commands = [self._addpath_commands()]
        # add the additional commands here
        if isinstance(matlab_args, BScanArguments):
            # handed BScanArgument, construct command from this
            self.matlab_commands += [matlab_args.create_bscan_argument()]
        else:
            # handed a list of strings, which are the individual matlab arguments
            self.matlab_commands += matlab_args
        # add the exit command because MATLAB does not give up control of the system after running a script >:(
        self.matlab_commands += [self._exit_command()]

    @staticmethod
    def _addpath_commands() -> str:
        """Produces a string whose text is a call to addpath() in matlab, so that matlab can find the run_bscan.m file and other supporting .m files needed for input generation."""
        command_str: str = "addpath("
        for extra_path in MATLAB_EXTRA_PATHS:
            command_str += "'" + extra_path + "',"
        # remove the extra comma after the last path was added
        command_str = command_str[:-1] + ")"
        return command_str

    @staticmethod
    def _exit_command() -> str:
        """Produces a string that corresponds to the exit command, which when run in matlab will close the instance."""
        return "exit"

    def create_command_string(self) -> str:
        """Sandwich each separate matlab command in the matlab_commands list into a single string matlab_command_string, which can be passed to matlab on the command line as:
        matlab -r matlab_command_string

        The string returned is matlab_command_string.
        """
        matlab_command_string = "; ".join(self.matlab_commands) + ";"
        return matlab_command_string

    def run(self) -> None:
        """Run the following matlab scripts/functions, in order, in THE SAME SESSION BECAUSE PATHING THINGS:

        addpath(EXTRA_MATLAB_PATHS)         <- run_bscan function & other required functions will be placed on the path
        run_bscan(<CORRECT CALL OPTIONS>)   <- populate call options with stuff from the config file
        exit                                <- prevent matlab from hanging after generating the data
        """
        # Generate the command-line call
        args = MATLAB_OPTS + [self.create_command_string()]
        # Run the matlab command
        p = Popen([MATLAB_PATH, *args], stdout=PIPE, cwd=LOCATION_OF_THIS_FILE)
        matlab_infodump, _ = p.communicate()
        # Save outputs for debugging later
        self.return_code = p.returncode
        self.stdout = matlab_infodump.decode()
        return


class GenerationData:
    """Handles the information and processes related to creating the input data for a particular test case.

    Information is read from a config.yaml file, and the contents of the "generation" field are validated against what is expected.

    From this information, the matlab_command is created. This command is the equivalent of running the old run_{pstd,fdtd}_bscan.m scripts on the test in question, however it is now sufficiently generalised so that we do not need an individual run file for each test.

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

    # The matlab command that generates the input data for this test
    matlab_command: MATLABCommand

    def __init__(self, config_file_location: Union[str, Path]) -> None:
        """Initialise instance by reading information from a config file"""
        # Fetch config file and options
        self._config_file_location = Path(config_file_location)
        config_file = open(config_file_location, "r")
        config = yaml.safe_load(config_file)
        config_file.close()

        # Get the ID of the test & thus the directory to save to
        self.test_id = config["test_id"]
        self.test_dir = Path(LOCATION_OF_THIS_FILE, ("arc_" + self.test_id))
        # Information about how this test's input data is generated
        self._generation_options = config["input_generation"]
        self._validate_config_options()

        # Setup the object based on the inputs
        self._setup_member_vars_from_config()

        # Setup the command from the member variables
        self.matlab_command = self._setup_matlab_command()
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

    def _setup_matlab_command(self) -> MATLABCommand:
        """Setup the matlab command that will generate the input data"""
        bscan_args = BScanArguments(self.test_dir, self.input_file)
        return MATLABCommand(bscan_args)

    def generate(self) -> int:
        """Generate the input data to the test, as specified by this instance's member values.

        Returns the exit code of the sequence of matlab commands that were run."""
        # ensure that the directory to place the output into exists, or create it otherwise
        self._find_or_create_test_dir()

        # generate the input data for this test
        self.matlab_command.run()

        # cleanup auxillary .mat files that are placed into this directory
        for aux_mat in sorted(glob(LOCATION_OF_THIS_FILE + "/*.mat")):
            os.remove(aux_mat)
        return self.matlab_command.return_code
