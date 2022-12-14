from pathlib import Path

import yaml
from utils import HDF5File, run_tdms

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

        # extract the information about this run
        # all required fields have been validated, so they exist
        for required_info in run_required_info:
            setattr(self, required_info, information[required_info])
        # attempt to extract optional fields, otherwise fill with default values
        for optional_info, default_option in run_optional_info:
            if optional_info in information.keys():
                # optional information provided
                setattr(self, optional_info, information[optional_info])
            else:
                # use default
                setattr(self, optional_info, default_option)
        return

    def get_output_file_name(self) -> None:
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

    def assert_output_matches_reference(self) -> None:
        # load reference file
        reference_file = HDF5File(self.reference)
        # load output file
        output_file = HDF5File(self.get_output_file_name())
        # assert file contents match
        assert output_file.matches(reference_file)
        return

    def run(self) -> None:
        tdms_command_call = self.create_tdms_call_options()
        # unpack into arguments to run_tdms
        print(f"Calling tdms with arguments {tdms_command_call}")
        run_tdms(*tdms_command_call)
        # assert the output matches the reference
        self.assert_output_matches_reference()
        return


class YAMLTestConfig:
    def __init__(self, config_file_path=Path("./config.yaml")) -> None:
        # get config information
        config = open_and_validate_config(config_file_path)

        # read test id
        self.test_id = config["test_id"]
        # read test information and create test runs
        self.generate_test_list(config["tests"])

        return

    def generate_test_list(self, runs: dict) -> list[TDMSSystemTest]:
        # check that all runs have the required information
        validate_run_information(runs)
        # unpack the test information
        self.run_list = []
        for run_name, run_info in runs.keys():
            self.run_list.append(TDMSSystemTest(run_name, run_info))

        return
