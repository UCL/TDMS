import os
from glob import glob
from pathlib import Path
from typing import Union

import yaml
from bscan_arguments import BScanArguments, MATLABEngineWrapper

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))


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
