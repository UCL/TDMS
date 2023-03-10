import os
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
from hdf5storage import loadmat, savemat

LOCATION_OF_THIS_FILE = os.path.abspath(os.path.dirname(__file__))


@dataclass
class OverwriteWith:
    """Small container class that will be used when overwriting the "flag" variables use_pstd and use_bli when a config file specifies an 'adjust' field.

    Intended workflow along the lines of:
    use_bli = OverwriteWith(whether_to_overwrite, value_to_replace_with)
    if use_bli.write:
        h5data['use_bli'] = use_bli.value
    else:
        pass
    """

    # Flags whether or not to overwrite a value with self.value
    write: bool
    # The value to write if we are intending to replace a value in a file
    value: Union[bool, int, float, None]


class MATFileOptionEdit:
    """Class that handles the generation of input .mat files from existing input .mat files. For certain system tests, we run tdms across the same numerical input data, but change the interpolation method and/or solver method flags. Since these are built into the input .mat files themselves, this requires a different input file for each run of tdms, despite the numerical data to be used being the same.

    To minimise the time spent regenerating the same numerical data, we can "adjust" the flag variables that control the interpolation/solver method by copying .mat input files and changing the value of the flag variables. We still produce multiple .mat inputs, but we only need one input_file.m and don't have to regenerate the numerical inputs each time.

    This class allows us to adjust the values of:
    - use_pstd: if true, we use pstd over fdtd
    - use_bli: if true, we use band-limited interpolation over cubic
    """

    # The test directory to search in and write to
    test_dir: str

    # The .mat file to produce
    to_produce: str
    # The .mat file to adjust
    adjust_file: str

    # Value of use_pstd to set to
    use_pstd: OverwriteWith
    # Value of use_bli to set to
    use_bli: OverwriteWith

    def __init__(
        self,
        test_directory: Union[Path, str],
        to_produce: Union[str, Path],
        options: dict[str, Union[str, bool, int]],
    ) -> None:
        """Create the object by passing in a field and its value(s) from the config file, assuming the options passed have the "adjust" field present."""
        self.test_dir = str(test_directory)

        # Get name of input file to produce immediately
        self.to_produce = str(Path(self.test_dir, to_produce))

        # Read the name of the .mat file to adjust
        self.adjust_file = str(Path(self.test_dir, options["adjust"]))

        # Read the variable values that we want to override, provided they exist
        # Check if we want to overwrite the solver method
        if "solver_method" in options.keys():
            # If this key is present, we want to overwrite the value of use_pstd with the value provided
            if options["solver_method"] == "pstd":
                # Overwrite use_pstd with 0 (pstd)
                self.use_pstd = OverwriteWith(True, True)
            elif options["solver_method"] == "fdtd":
                # Overwrite use_pstd with 1 (fdtd, default)
                self.use_pstd = OverwriteWith(True, False)
            else:
                raise RuntimeError(
                    f"Error: {options['solver_method']} is not a valid solver method"
                )
        else:
            # This key is not present, so we do not want to overwrite
            # Note that the self.use_pstd.value variable will not be used in this case, so we can set it to anything
            self.use_pstd = OverwriteWith(False, None)

        # Check if we want to overwrite the interpolation method
        if "interpolation" in options.keys():
            # We will be overwriting the value here
            if options["interpolation"] == "bli":
                # Overwrite use_bli with 2 (band-limited)
                self.use_bli = OverwriteWith(True, True)
            elif options["interpolation"] == "cubic":
                # Overwrite use_bli with 1 (cubic, default)
                self.use_bli = OverwriteWith(True, False)
            else:
                raise RuntimeError(
                    f"Error: {options['interpolation']} is not a valid interpolation method"
                )
        else:
            # Overwrite not requested
            self.use_bli = OverwriteWith(False, None)

        # Setup complete, return
        return

    def adjust(self) -> None:
        """Perform the adjustments to the .mat input file, and produce the output file requested."""
        # Load the .mat file that is to be adjusted
        mat_data = loadmat(self.adjust_file)
        # Adjust use_pstd if necessary
        if self.use_pstd.write:
            mat_data["use_pstd"] = np.array([[self.use_pstd.value]], dtype=bool)
        # Adjust use_bli if necessary
        if self.use_bli.write:
            mat_data["use_bli"] = np.array([[self.use_bli.value]], dtype=bool)
        # Write the data to the desired location - THROWS DEPRECATION WARNING FROM NUMPY!
        savemat(self.to_produce, mat_data, format="7.3", store_python_metadata=False)
        # Complete, return
        return
