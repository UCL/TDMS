import os
from typing import Union

import matlab.engine as matlab
from bscan_arguments import BScanArguments
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


class MATLABEngineWrapper:
    """When we regenerate input data, we always need to add the bscan/ and matlab/ directories to the MATLAB instance's search path. correspondingly, we always want to kill the MATLAB instance after running the run_bscan function and generating the data.

    This class serves as a wrapper for that purpose. It stores instance(s) of the BScanArguments class, which it will run in sequence between the aforementioned addpath() setup and then engine shutdown. The .run() method performs exactly this.
    """

    # List containing the subsequent matlab commands to be executed by the interpreter
    bscan_calls: list[BScanArguments]

    # The MATLAB session that will run
    engine: MatlabEngine
    # The working directory this engine operated in. Useful for cleanup operations later
    cwd: str

    def __init__(
        self,
        bscans: Union[BScanArguments, list[BScanArguments]],
        engine_cwd=LOCATION_OF_THIS_FILE,
    ) -> None:
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
        # record the working directory for the engine to be started in, though
        self.cwd = engine_cwd
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
            # Move to the requested working directory
            self.engine.cd(self.cwd)
            # Add necessary paths
            self._addpath_commands()

            # Run every bscan call
            for b in self.bscan_calls:
                b.run_bscan(self.engine)

            # If requested, kill the instance
            if kill_on_complete:
                self._stop_engine()
        return
