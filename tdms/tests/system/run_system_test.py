import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from warnings import warn
from zipfile import ZipFile

import yaml
from pytest_check import check
from utils import HDF5File, Result, run_tdms

LOCATION_OF_THIS_FILE = Path(os.path.abspath(os.path.dirname(__file__)))
ZIP_DIR = LOCATION_OF_THIS_FILE / "data"


@dataclass
class TDMSRun:
    """One run, or execution, of the TDMS executable. That is, one call to

    tdms [OPTIONS] [input_file] [gridfile] [output_file]

    that is required as part of a single system test.
    The execute() command executes the above call to TDMS
    """

    # ID of this run for diagnostics/debug
    id: str

    # The .mat input to the tdms executable
    input_file: Path
    # Output .mat file to write to
    output_file: Path
    # Gridfile .mat input, if provided
    gridfile: Path | None
    # Additional command-line flags
    flags: list[str]
    # Command in list form to be fed into utils.run_tdms
    command: list[str] = field(init=False)

    # The reference data file (relative to the top-level of the .zip archive for this test)
    ref_data: Path

    def __post_init__(self):
        """Setup the command that will be passed to utils.run_tdms."""
        self.command = self.flags + [str(self.input_file)]
        # Gridfile comes next if it was provided
        if self.gridfile:
            self.command.append(str(self.gridfile))
        # Output file is final input on command line
        self.command.append(str(self.output_file))

    def __str__(self) -> str:
        return f"{' '.join(['tdms'] + self.command)}"

    def execute(self) -> Result:
        """Call tdms with the commands specified in this instance"""
        sys.stdout.write(f"Now running:\n\t {self}")
        start = time.time()
        result = run_tdms(*self.command)
        ellapsed_s = time.time() - start
        sys.stdout.write(f"\t (took {ellapsed_s}s)")
        return result


def run_system_test(config_filepath: Path | str) -> dict[str, bool]:
    """Perform the system test detailed in the config file passed.

    :param config_filepath: Path to a .yaml file containing the information about the system test to run.
    :returns: Dictionary indexed by the run_ids. Values are boolean; True if that run passed, and False if the run failed.
    """
    with open(config_filepath, "r") as file:
        test_information = yaml.safe_load(file)

    # Track the test's id (arc_{test_id}) and the information about the runs it contains
    test_id = test_information["test_id"]

    # Infer the location of the input data and .zip archive
    input_data_folder = (
        LOCATION_OF_THIS_FILE / "data" / "input_generation" / f"arc_{test_id}"
    )
    ref_data_zip_foler = ZipFile(ZIP_DIR / f"arc_{test_id}.zip", "r")

    # Infer the names of the input files that make up this run
    mat_inputs = [mfile for mfile in test_information.keys() if mfile != "test_id"]

    # Setup the runs that are part of this test, using the test_dict
    tdms_runs: list[TDMSRun] = []
    # These are reference files that have been pulled out of .zip archives and will need cleaning up afterwards
    # set() ensures that duplicate filenames will not be added if already present
    extracted_ref_files = set()
    for mat_input_file in mat_inputs:
        # If this input file is not used in any runs (it always should be, but it never hurts to check), then continue to the next element without running any tests
        if "runs" in test_information[mat_input_file].keys():
            # Perform each run that uses this file as it's input
            path_to_input_m_file = input_data_folder / (
                os.path.splitext(mat_input_file)[0] + ".mat"
            )

            runs = test_information[mat_input_file]["runs"]
            # The dict itself contains many runs, so we need to loop again
            for run_name, run_info in runs.items():
                # Create a suitable name for the output of the run, and recycle the input_data_folder for the location of the output
                output_file = input_data_folder / f"output_{test_id}_{run_name}.mat"

                # Fetch the reference output that this run needs to compare to
                ref_output = run_info["reference"]
                extracted_ref_files.add(ref_output)

                # Fetch the gridfile, if it exists
                if "gridfile" in run_info.keys():
                    gridfile = run_info["gridfile"]
                else:
                    gridfile = None

                # Determine if there are any flags that need to be passed to the tdms run
                run_flags = []
                # For each flag that we want to pick up on from the config file,
                # the following syntax can be used.
                # if "config_member_indicating_flag" in run_info and run_info["config_member_indicating_flag"]:
                #     run_flags += ["CLI_flag_string"]

                # Create the run and store it with the reference data file that it needs
                tdms_runs.append(
                    TDMSRun(
                        run_name,
                        path_to_input_m_file,
                        output_file,
                        gridfile,
                        run_flags,
                        ref_output,
                    )
                )
        else:
            warn(f"Input file {mat_input_file} corresponds to no TDMS runs")

    # Extract all unique reference data files from the .zip archive
    dir_to_extract_refs_to = str(input_data_folder / "reference_data")
    for ref_data in extracted_ref_files:
        # extract() cannot take Path arguments
        ref_data_zip_foler.extract(ref_data, dir_to_extract_refs_to)

    # Perform each run and log the success, assume runs fail by defualt
    run_success = dict()

    # Perform each run, and compare the generated output to the reference output
    for run in tdms_runs:
        # Run the TDMS call
        run.execute()

        # Compare the output to the reference data
        run_output = HDF5File(run.output_file)
        ref_output = HDF5File(f"{dir_to_extract_refs_to}/{run.ref_data}")

        # Use check, so that all runs are performed before we report test failure
        # This also ensures that we cleanup the reference data REGARDLESS of test success
        with check:
            test_passed, comparison_information = run_output.matches(
                ref_output, return_message=True
            )
            assert test_passed, f"In {test_id} -> {run.id}: {comparison_information}"
        # Append whether or not the test passed
        run_success[run.id] = test_passed

    # Cleanup the .mat files we had to copy across
    for mat_file in extracted_ref_files:
        os.remove(f"{dir_to_extract_refs_to}/{mat_file}")
    # then remove the (what should be empty) directory
    os.rmdir(dir_to_extract_refs_to)

    return run_success
