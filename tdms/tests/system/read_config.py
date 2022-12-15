from pathlib import Path

import yaml
from utils import run_tdms

# config files must provide this information
required_fields = ["test_id", "tests"]
# any runs must provide this information
run_required_info = ["input_file", "reference"]
# runs may also specify this information, but if they don't we use these defaults
run_optional_info = {"fdtd_solver": False, "cubic_interpolation": False}


def open_and_validate_config(config_file_path) -> dict:
    # open the config file
    config_file = open(config_file_path, "r")
    # read data into dictionary
    config = yaml.safe_load(config_file)
    # can safely close config file after reading data
    config_file.close()

    for required_field in required_fields:
        if required_field not in config.keys():
            raise RuntimeError(f"{required_field} not provided in {config_file_path}")

    # pass the data from the file back
    return config


def validate_run_information(runs: dict) -> None:
    for test_id, test_info in runs.items():
        for required_info in run_required_info:
            if required_info not in test_info.keys():
                raise RuntimeError(f"{required_info} not provided in {test_id}")
    return


class TDMSSystemTest:
    def __init__(self, name: str, information: dict) -> None:
        # associate the command to the run name
        self.run_id = name

        # check that the user has not provided any unexpected information
        for info in information.keys():
            if (info not in run_required_info) and (
                info not in run_optional_info.keys()
            ):
                # we don't expect this information! Throw an error
                raise RuntimeError(f"{self.run_id} got unexpected information: {info}")

        # extract the information about this run
        # all required fields have been validated, so they exist
        for required_info in run_required_info:
            setattr(self, required_info, information[required_info])
        # attempt to extract optional fields, otherwise fill with default values
        for optional_info, default_option in run_optional_info.items():
            if optional_info in information.keys():
                # optional information provided
                setattr(self, optional_info, information[optional_info])
            else:
                # use default
                setattr(self, optional_info, default_option)
        return

    def get_reference_file_name(self) -> str:
        return self.reference

    def get_output_file_name(self) -> str:
        # output files generated from running tdms are named as below
        # these files are to be compared to self.reference
        return f"./{self.run_id}_output.mat"

    def create_tdms_call_options(self) -> list[str]:
        tdms_call_options = []
        # add optional flags
        if self.fdtd_solver:
            tdms_call_options.append("--finite-difference")
        if self.cubic_interpolation:
            tdms_call_options.append("-c")

        # add input/output file call TODO: what if we want to specify a gridfile (no tests currently do this)
        tdms_call_options.append(self.input_file)
        # add output file to tdms call
        tdms_call_options.append(self.get_output_file_name())

        # return a list that can be passed to run_tdms with * expansion
        return tdms_call_options

    def run(self) -> None:
        tdms_command_call = self.create_tdms_call_options()
        # unpack into arguments to run_tdms
        print(f"Calling tdms with arguments {tdms_command_call}")
        run_tdms(*tdms_command_call)
        return


class YAMLTestConfig:
    """
    Each test has a config file, config.yaml, at the top-level of the .zip folder.
    This config file contains information about the tests that are to be run using the data and reference outputs contained in the folder.

    Config file syntax is presently:
    test_id: XX
    tests:
        test_name_no_spaces:
            input_file: relative_path_to_input
            reference: relative_path_to_reference_output
            fdtd_solver: True/False [optional, default is False is absent]
            cubic_interpolation: True/False [optional, default is False if absent]

    You can define multiple tests using the usual yaml syntax. For example, arc_01.zip contains 4 tests which are a mixture of cubic and band-limited interpolation tests.

    test_id: 01
    tests:
        fs_bli:
            input_file: pstd_fs_input.mat
            reference: pstd_fs_bli_reference.mat
        fs_cubic:
            input_file: pstd_fs_input.mat
            reference: pstd_fs_cubic_reference.mat
            cubic_interpolation: True
        cyl_bli:
            input_file: pstd_cyl_input.mat
            reference: pstd_cyl_bli_reference.mat
        cycl_cubic:
            input_file: pstd_cyl_input.mat
            reference: pstd_cyl_bli_reference.mat
            cubic_interpolation: True

    It is also possible to mix-and-match the fdtd_solver in a similar way.
    """

    def __init__(self, config_file_path=Path("./config.yaml")) -> None:
        # get config information
        config = open_and_validate_config(config_file_path)

        # read test id
        self.test_id = config["test_id"]
        # read test information and create test runs
        self.run_list = self.generate_test_list(config["tests"])

        return

    def generate_test_list(self, runs: dict) -> list[TDMSSystemTest]:
        # check that all runs have the required information
        validate_run_information(runs)
        # unpack the test information
        run_list = []
        for run_name, run_info in runs.items():
            run_list.append(TDMSSystemTest(run_name, run_info))

        return run_list
