import os

import matlab.engine

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))
PY_INPUT_REGEN_LOCATION = os.path.abspath(LOCATION_OF_THIS_FILE + "../")

# Additional options for running matlab on the command-line
MATLAB_OPTS_LIST = ["-nodisplay", "-nodesktop", "-nosplash", "-r"]
MATLAB_STARTUP_OPTS = " ".join(MATLAB_OPTS_LIST)

# Paths to matlab functions not in LOCATION_OF_THIS_FILE
MATLAB_EXTRA_PATHS = [
    os.path.abspath(LOCATION_OF_THIS_FILE + "../../../data/input_generation/matlab"),
    os.path.abspath(LOCATION_OF_THIS_FILE + "../../../data/input_generation/bscan"),
]


def create_engine_for_testing(
    working_dir=PY_INPUT_REGEN_LOCATION,
) -> matlab.engine.MatlabEngine:
    """Creates a MatlabEngine instance runing in the directory provided (default is the py_input_regen directory), and with the following directories added to the MATLABPATH:
    - tdms/tests/system/data/input_generation/matlab
    - tdms/tests/system/data/input_generation/bscan

    The user is responsible for quit()-ing the engine instance.

    :param working_dir: The working directory that the engine should be launched in.
    :returns: An active MatlabEngine instance with the setup described.
    """
    engine = matlab.engine.start_matlab(MATLAB_STARTUP_OPTS)
    # Append additional paths
    engine.cd(working_dir)
    for PATH in MATLAB_EXTRA_PATHS:
        engine.addpath(PATH)

    return engine
