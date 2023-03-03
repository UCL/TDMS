from dataclasses import dataclass
from pathlib import Path
from typing import Union

from utils import run_tdms


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
