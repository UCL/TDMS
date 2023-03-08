import os
from pathlib import Path
from typing import Union
from warnings import warn
from zipfile import ZipFile

import yaml
from pytest_check import check
from tdms_run import TDMSRun, TDMSRunAndReference
from utils import HDF5File

LOCATION_OF_THIS_FILE = Path(os.path.abspath(os.path.dirname(__file__)))
ZIP_DIR = LOCATION_OF_THIS_FILE / "data"


class TDMSSystemTest:
    """Instance of one tdms system test. This object controls all the tdms executations that required for the system test arc_{test_id} to run, and handles the setup and tear-down of each run.

    The runs themselves are performed by TDMSRun instances.
    """

    # test_id of this test. This determines filenames and other important locations
    test_id: str
    # Folder in which locally generated input data will be placed. We will also recycle this directory as where the output of the tdms runs, and the reference data, will be copied to before cleanup
    input_data_folder: Path
    # The .zip file that contains reference data
    ref_data_zip_foler: ZipFile

    # A list of all tdms runs that are part of this system test, and their asociated reference files to X-reference against
    tdms_runs: list[TDMSRunAndReference]
    # List of the reference data files we will be pulling out of the .zip archive (and cleaning up after)
    ref_files: list[str]

    def __init__(self, config_file: Union[Path, str]) -> None:
        """Initialise the system test by passing in the confif_{test_id}.yaml file."""
        # Load data from config file
        with open(config_file, "r") as opened_config_file:
            # Read data into dictionary, and take the "tests" value
            config = yaml.safe_load(opened_config_file)

        # Self-assign test_id field
        self.test_id = config["test_id"]
        # From this, infer the location of the input data and .zip archive
        self.input_data_folder = (
            LOCATION_OF_THIS_FILE / "data" / "input_generation" / f"arc_{self.test_id}"
        )
        self.ref_data_zip_foler = ZipFile(ZIP_DIR / f"arc_{self.test_id}.zip", "r")

        # Setup the runs that are part of this test, using the test_dict.
        # Existence of files, etc is checked by the TDMSRun object.
        # The non-"test_id" keys in config are names of input files, which contain a "runs" key that details the runs they are used in.
        input_files = list(config.keys())
        input_files.remove("test_id")
        self.ref_files = []
        self.tdms_runs = []
        for mat_input_file in input_files:
            # If this input file is used in any runs (it always should, but it never hurts to check)
            if "runs" in config[mat_input_file].keys():
                # Setup the runs that use this input file
                self._setup_runs_that_use_input(
                    mat_input_file, config[mat_input_file]["runs"]
                )
            else:
                warn(f"Input file {mat_input_file} corresponds to no TDMS runs")

        # Setup is now complete, return
        return

    def _setup_runs_that_use_input(
        self, mat_input: str, runs: dict[str, dict[str, str]]
    ) -> None:
        """Add the runs that make use of the given .mat file (mat_input) to the list of runs that need to be executed.

        runs should be the value of the config_{test_id}.yaml[mat_input]["runs"] field.
        """
        # For safety with absolute paths, ensure that the .mat extension is present on the input file name
        # The input file for each of these runs is otherwise the same between runs
        input_file = self.input_data_folder / Path(
            os.path.splitext(mat_input)[0] + ".mat"
        )

        # The dict itself contains many runs, so we need to loop again
        for run_name, run_info in runs.items():
            # Fetch the reference output that this run needs to compare to
            ref_output = run_info["reference"]
            if ref_output not in self.ref_files:
                # This is a new reference file that will be pulled out of the .zip archive, and which needs cleaning up afterwards
                self.ref_files.append(ref_output)

            # Create a suitable name for the output of the run, and recycle the input_data_folder for the location of the output
            output_file = self.input_data_folder / self._choose_output_name(run_name)

            # Fetch the gridfile, if it exists
            if "gridfile" in run_info.keys():
                gridfile = run_info["gridfile"]
            else:
                gridfile = None

            # Determine if there are any flags that need to be passed to the tdms run
            run_flags = []
            if "cubic_interpolation" in run_info.keys():
                if run_info["cubic_interpolation"]:
                    run_flags += ["-c"]
            if "fdtd_solver" in run_info.keys():
                if run_info["fdtd_solver"]:
                    run_flags += ["--finite-difference"]

            # Create the run and store it with the reference data file that it needs
            self.tdms_runs.append(
                TDMSRunAndReference(
                    TDMSRun(run_name, input_file, output_file, gridfile, run_flags),
                    ref_output,
                )
            )
        # We have now appended all the runs that use mat_input as the input file to the list of runs to perform
        return

    def _choose_output_name(self, run_name: str) -> str:
        """Creates a name for the (temporary) output file from a call to tdms, given the {run_name}."""
        return f"output_{self.test_id}_{run_name}.mat"

    def _dir_to_extract_refs_to(self) -> str:
        """Provides the name of the directory into which we will be extracting reference data from the .zip archives, to compare to the output of the latest test runs."""
        return str(self.input_data_folder / "reference_data")

    def _ref_file_cleanup(self) -> None:
        """Remove any .mat output files that we had to unzip for comparison with the test runs"""
        for mat_file in self.ref_files:
            os.remove(f"{self._dir_to_extract_refs_to()}/{mat_file}")
        # then remove the (what should be empty) directory
        os.rmdir(self._dir_to_extract_refs_to())
        return

    def perform_all_runs(self, cleanup_ref_outputs=True) -> list[bool]:
        """Execute all calls to tdms mandated by the runs in this system test, and compare the results to the reference data.

        Return a list of bools: True entries indicate a passing run, False entries indicate a failing run
        """
        # Extract the reference data from the .zip archive
        for ref_data in self.ref_files:
            self.ref_data_zip_foler.extract(ref_data, self._dir_to_extract_refs_to())

        # The output - a list indicating which runs passed / failed
        run_passes = []

        # Perform each run, and compare the generated output to the reference output
        for run in self.tdms_runs:
            # Run the TDMS call
            exit_code = run.tdms_run.run()

            # Compare the output to the reference data
            run_output = HDF5File(run.tdms_run.output_file)
            ref_output = HDF5File(f"{self._dir_to_extract_refs_to()}/{run.ref_data}")
            # Use check, so that all runs are performed before we report test failure
            # This also ensures that we cleanup the reference data REGARDLESS of test success
            with check:
                test_passed, comparison_information = run_output.matches(
                    ref_output, return_message=True
                )
                assert (
                    test_passed
                ), f"In {self.test_id} -> {run.tdms_run.run_id}: {comparison_information}"
            # Append whether or not the test passed
            run_passes.append(test_passed)

        # Cleanup the .mat files we had to copy across
        if cleanup_ref_outputs:
            self._ref_file_cleanup()

        return run_passes
