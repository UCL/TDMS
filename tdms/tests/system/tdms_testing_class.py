import os
from dataclasses import dataclass
from pathlib import Path
from typing import Union
from zipfile import ZipFile

from pytest_check import check
from utils import HDF5File, run_tdms

LOCATION_OF_THIS_FILE = Path(os.path.abspath(os.path.dirname(__file__)))
ZIP_DIR = LOCATION_OF_THIS_FILE / "data"


class TDMSRun:
    """One run, or execution, of the TDMS executable. That is, one call to

    tdms [OPTIONS] [input_file] [gridfile] [output_file]

    that is required as part of a single system test.
    The run() command executes the above call to TDMS
    """

    # The ID of this run within the wider system test
    run_id: str

    # The .mat input to the tdms executable
    input_file: str
    # Output .mat file to write to
    output_file: str
    # Gridfile .mat input, if provided
    gridfile: str
    # Additional command-line flags
    flags: list[str]

    # stdout produced by the call to tdms
    run_stdout: str

    def __init__(
        self,
        name: str,
        input_file: Union[Path, str],
        output_file: Union[Path, str],
        gridfile: Union[Path, str, None] = None,
        flags: list[str] = [],
    ):
        """Initialise having been passed a member of config_{test_id}.yaml[tests].
        That is, run_information = config_{test_id}.yaml[tests][name].
        """
        # Required: logging and debugging
        self.run_id = name

        # Required: executation of tdms
        # Input file
        self.input_file = str(input_file)
        if not Path(self.input_file).exists():
            raise RuntimeError(f"Input file {self.input_file} not found.")
        # Location of output file
        self.output_file = str(output_file)

        # Optional: executation of tdms
        # Gridfile, if it is passed
        if gridfile:
            self.gridfile = str(gridfile)
            if not Path(self.gridfile).exists():
                raise RuntimeError(f"Gridfile {self.gridfile} not found.")
        else:
            self.gridfile = ""
        # Command-line flags
        self.flags = flags
        return

    def __str__(self) -> str:
        return f"{' '.join(['tdms'] + self._assemble_command())}"

    def __repr__(self) -> str:
        return f"Instance of TDMSRun holding the command: \n\t{' '.join(['tdms'] + self._assemble_command())}"

    def _assemble_command(self) -> list[str]:
        """Generates a list of strings which form the command-line input to tdms, for this run."""
        # List of command-line arguments that are to be passed to tdms executable
        # Always starts with the flags and the input file
        tdms_cl_inputs = self.flags + [self.input_file]
        # Gridfile comes next if it was provided
        if self.gridfile:
            tdms_cl_inputs.append(self.gridfile)
        # Output file is final input on command line
        tdms_cl_inputs.append(self.output_file)
        # Return list of commands
        return tdms_cl_inputs

    def run(self, print_command: bool = True) -> int:
        """Run tdms according to the specifications of this instance. print_command controls whether or not to log the command actually being executed to stdout.

        Return the exit code from the shell subprocess that was spawned. stdout is placed into the run_stdout member.
        """
        # List of command-line arguments that are to be passed to tdms executable
        tdms_cl_inputs = self.flags + [self.input_file]
        # Gridfile comes next if it was provided
        if self.gridfile:
            tdms_cl_inputs.append(self.gridfile)
        # Output file is final input on command line
        tdms_cl_inputs.append(self.output_file)

        # Print command before run if requested
        if print_command:
            print("Now running:\n\t", self)
        # Run executable and report return code
        tdms_result = run_tdms(*tdms_cl_inputs)
        # Preserve stdout in case we need to debug
        self.run_stdout = tdms_result.stdout
        # Return the exit code of the shell subprocess
        return tdms_result.return_code


@dataclass
class TDMSRunAndReference:
    tdms_run: TDMSRun  # The object that handles the tdms call itself
    ref_data: str  # The name of the reference data within the .zip archive


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

    def __init__(self, test_id: str, run_dict: dict[str, dict[str, str]]) -> None:
        """Initialise the system test by passing in the information from config_{test_id}.yaml["tests"]."""
        self.test_id = test_id
        # From this, infer the location of the input data and .zip archive
        self.input_data_folder = (
            LOCATION_OF_THIS_FILE / "data" / "input_generation" / f"arc_{self.test_id}"
        )
        self.ref_data_zip_foler = ZipFile(ZIP_DIR / f"arc_{self.test_id}.zip", "r")

        # Setup the runs that are part of this test, using the test_dict
        # Existence of files, etc is checked by the TDMSRun object
        self.ref_files = []
        self.tdms_runs = []
        for run_name, run_info in run_dict.items():
            # Fetch the input file
            input_file = self.input_data_folder / run_info["input_file"]

            # Create a suitable name for the output of the run, and recycle the input_data_folder for the location of the output
            output_file = self.input_data_folder / self._choose_output_name(run_name)
            # Fetch the filename of the reference data that this run compares it's output to
            ref_output = run_info["reference"]
            if ref_output not in self.ref_files:
                # This is a new reference file that will be pulled out of the .zip archive, and which needs cleaning up afterwards
                self.ref_files.append(ref_output)

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

    def perform_all_runs(self) -> list[bool]:
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
        self._ref_file_cleanup()

        return run_passes
