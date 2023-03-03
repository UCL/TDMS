import os
from glob import glob
from pathlib import Path
from typing import Union

import yaml
from bscan_arguments import BScanArguments
from matlab_engine_wrapper import MATLABEngineWrapper

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

    # ID of the test we are generating input data for
    test_id: str
    # Absolute path to the directory into which the input data should be placed
    test_dir: Path

    # The matlab instance that generates the input data for this test
    matlab_instance: MATLABEngineWrapper
    # The .mat input files that we will produce
    matfiles_to_produce: list[str]

    def __init__(
        self,
        config_file_location: Union[str, Path],
        test_dir: Union[Path, str] = LOCATION_OF_THIS_FILE,
    ) -> None:
        """Initialise instance by reading information from a config file"""
        # Fetch config file and options
        self._config_file_location = Path(config_file_location)
        with open(config_file_location, "r") as config_file:
            config = yaml.safe_load(config_file)

        # Get the ID of the test & thus the directory to save to
        self.test_id = config["test_id"]
        self.test_dir = Path(test_dir, f"arc_{self.test_id}")

        # Determine the .mat files that we will be generating
        # Non-"test_id" fields are the .mat input file names
        self.matfiles_to_produce = list(config.keys())
        self.matfiles_to_produce.remove("test_id")

        # Setup the command from the member variables
        self.matlab_instance = self._setup_matlab_instance(config)
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

    def _setup_matlab_instance(self, config: dict[str, any]) -> MATLABEngineWrapper:
        """Setup the BScan commands that will be run by the engine, using the information in the config file."""
        # This is the list of BScan commands that this config file wants us to run
        bscan_list = []
        for mat_file in self.matfiles_to_produce:
            bscan_list.append(BScanArguments(self.test_dir, mat_file, config[mat_file]))
        # Return a matlab engine that is ready to run each of these commands
        return MATLABEngineWrapper(bscan_list, str(self.test_dir))

    def _cleanup(self) -> None:
        """Cleanup auxillary .mat files that are placed into the working directory by the run_bscan function."""
        # pull the working directory of the matlab engine for cleanup reasons
        matlab_cwd = self.matlab_instance.cwd
        # Create set of all files to cleanup
        mat_files_dumped_here = set(glob(str(self.test_dir) + "/*.mat"))
        mat_files_dumped_cwd = set(glob(matlab_cwd + "/*.mat"))
        dumped_mat_files = mat_files_dumped_cwd | mat_files_dumped_here
        # exclude the input files themselves from deletion
        matfiles = [
            f"{self.test_dir}/{matfile}.mat" for matfile in self.matfiles_to_produce
        ]
        dumped_mat_files = dumped_mat_files - set(matfiles)
        # purge .mat files
        for aux_mat in dumped_mat_files:
            os.remove(aux_mat)
        return

    def generate(self) -> None:
        """Generate the input data to the test, as specified by this instance's member values.

        Returns the exit code of the sequence of matlab commands that were run."""
        # ensure that the directory to place the output into exists, or create it otherwise
        self._find_or_create_test_dir()

        # generate the input data for this test
        self.matlab_instance.run()

        # cleanup auxillary .mat files
        self._cleanup()
        return
