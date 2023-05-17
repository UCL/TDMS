#########
#
# See YAMLTestConfig.help() for details on writing config files for tests.
#
#########
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
    """Open the config file at config_file_path and validate the specifications it provides.

    Config files are required to have all required_fields present.
    """
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


def validate_run_information(runs: dict[str, dict[str, str]]) -> None:
    """Validates the run (call to tdms) specifies all required fields.

    Each entry of "runs" contains the specifications for one run of the tdms executable.
    Each such "run" is a dict that is required to have a key corresponding to the entries of run_required_info.
    run_optional_info contains other properties that might be specified by a key in "run", and the defaults to use if these properties are not provided in the config file.
    """
    for test_id, test_info in runs.items():
        for required_info in run_required_info:
            if required_info not in test_info.keys():
                raise RuntimeError(f"{required_info} not provided in {test_id}")
    return


class TDMSSystemTest:
    """Objects of this class handle a single "run" specified in a config file, producing the command-line call to tdms and then running said command."""

    def __init__(self, name: str, information: dict) -> None:
        """Given the name of a run (EG: fs_bli) and the corresponding information about it from the config file, setup the instance."""
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
        # all required fields have been validated above, so they exist
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
        """Return the name of the reference file to compare the output of this call to tdms to."""
        # reference is a required field, and an error is thrown if it is not specified. So this attribute always exists.
        return self.reference

    def get_output_file_name(self) -> str:
        """Returns the path that the output of the tdms executable will be saved to - this will be compared to self.reference."""
        # output files generated from running tdms are named as below
        # these files are to be compared to self.reference
        return f"./{self.run_id}_output.mat"

    def create_tdms_call_options(self) -> list[str]:
        """Returns a list of command-line options (in appropriate order) that forms the command-line call to tdms that this instance encodes.

        This output can be passed to utils.run_tdms to run the tdms executable with the appropriate options.
        """
        tdms_call_options = []
        # add optional CLI flags that need to be passed to TDMS here, by appending to tdms_call_options

        # add input/output file call TODO: what if we want to specify a gridfile (no tests currently do this)
        tdms_call_options.append(self.input_file)
        # add output file to tdms call
        tdms_call_options.append(self.get_output_file_name())

        # return a list that can be passed to run_tdms with * expansion
        return tdms_call_options

    def run(self) -> None:
        """Run the tdms command specified by this instance."""
        # create the options for the tdms call
        tdms_command_call = self.create_tdms_call_options()
        # display tdms call to the user for debugging/logging purposes
        print(f"Calling tdms with arguments {tdms_command_call}")
        # call run_tdms with these options
        run_tdms(*tdms_command_call)
        return


class YAMLTestConfig:
    """Parser class for config.yaml files contained in .zip folders containing system test data. See the help() method for more information."""

    def __init__(self, config_file_path=Path("./config.yaml")) -> None:
        """Parse the configuration file at config_file_path."""
        # get config information
        config = open_and_validate_config(config_file_path)

        # read test id
        self.test_id = config["test_id"]
        # read test information and create test runs
        self.run_list = self.generate_test_list(config["tests"])

        return

    def generate_test_list(self, runs: dict) -> list[TDMSSystemTest]:
        """Generate a list of system tests requested in the configuration file.

        A single .zip folder, and thus configuration file, specifies multiple system tests, or "runs" of the tdms executable. Having parsed the configuration file, we can now generate a list of TDMSSystemTest objects that corresponds to these "runs", with one TDMSSystemTest object per run specified.
        """
        # check that all runs have the required information
        validate_run_information(runs)
        # unpack the test information
        run_list = []
        for run_name, run_info in runs.items():
            run_list.append(TDMSSystemTest(run_name, run_info))
        # return the list of runs (system tests)
        return run_list

    def help(self) -> None:
        """Prints details about input file specification to the screen."""
        help_message = """
        A single .zip file of system test data may contain multiple system tests, referred to as "runs" in this codebase so that pytest doesn't try to go crazy. Each "run" consists of a single call to the tdms executable, using a particular input file and command-line options, and then a comparison to some reference data.

        So that pytest knows which input files and options match up to particular outputs, and the tdms options that are necessary to generate them, each .zip folder contains a config.yaml with this information. These config files are stored at the top-level of the .zip folder, and use the name "config.yaml" (although there is functionality to handle using a different name if someone so wishes).

        Config file syntax is:
        test_id: XX
        tests:
            test_name_1:
                input_file: input_path_1
                reference: output_path_1
                fdtd_solver: True/False
                cubic_interpolation: True/False
            test_name_2:
                input_file: input_path_2
                reference: output_path_1
                fdtd_solver: True/False
                cubic_interpolation: True/False
        Indentation is used to define blocks associated to the attribute one-indentation level higher.

        test_id is a string (although integers can also be supplied and will be converted to strings) that provides a unqiue identifier for the .zip folder. The zip folders typically follow the naming convention arc_{test_id}.zip, EG arc_01.zip and arc_example_fdtd.zip.

        Members under the "tests" identifier correspond to individual runs.

        The "test_id" and "tests" fields are required, failing to provide these will result in an error being thrown.

        Each run (member of "tests") is identified by its test_name, which must not contain spaces. These typically follow a naming convention of {solver_method}_{spatial_object} (or the reverse); for example pstd_fs or cyl_fdtd. The name can be anything you like (it is not actually used by the interpreter other than for debugging readability purposes) but for this reason we advise against naming all the runs "foobar".

        Within each run are the fields input_file, reference, and (optionally) cubic_interpolation and fdtd_solver.
        The latter two are bools indicating whether the corresponding tdms option (use cubic interpolation or the fdtd solver respectively) should be toggled on. If one of these fields is missing, the default value is False (use bandlimited interpolation or pstd solver methods respectively).
        The input_file and reference fields are compulsory, and must specifiy the path (RELATIVE to the top level of the zipped directory) to the input file to be passed to the tdms executable, and the reference output to compare the generated output to. The generated output is named automatically.

        An example config file is below:

        '''yaml
        test_id: config_example
        tests:
            pstd_and_bli:
                input_file: input_a
                reference: reference_a
            pstd_and_cubic:
                input_file: input_a
                reference: reference_b
                cubic_interpolation: True
            fdtd_and_bli:
                input_file: input_b
                reference: reference_b
                fdtd_solver: True
            fdtd_and_cubic:
                input_file: input_c
                reference: reference_c
                cubic_interpolation: True
                fdtd_solver: True
        '''

        This config file corresponds to 4 runs, the tdms commands being:
        - tdms input_a automatically_generated_output_name_1
        - tdms -c input_a automatically_generated_output_name_2
        - tdms --finite-difference input_b automatically_generated_output_name_3
        - tdms -c --finite-difference input_c automatically_generated_output_name_4

        And then the comparisons made being:
        - automatically_generated_output_name_1 VS reference_a
        - automatically_generated_output_name_2 VS reference_b
        - automatically_generated_output_name_3 VS reference_b
        - automatically_generated_output_name_4 VS reference_c
        """
        print(help_message)
        return
