from utils import HDF5File, download_data, run_tdms, work_in_zipped_dir
from zip_file_config_object import ZipFileConfig

# Config options for this test
id = 1
zenodo_link = "https://zenodo.org/record/6838866/files/arc_01.zip"
other_spatial_object = "cyl"
uses_fdtd = False

# Create the config object
test_config = ZipFileConfig(id, zenodo_link, other_spatial_object, uses_fdtd)
# Fetch data if it is not already available
if not test_config.zip_path.exists():
    download_data(test_config.zenodo_link, to=test_config.zip_path)


@work_in_zipped_dir(test_config.zip_path)
def test_inputs_in_zip() -> None:
    """
    Tests both cubic and band-limited interpolation outputs computed from the relevent input files in this .zip folder, against the corresponding reference outputs in this .zip folder.
    """

    for test_case in test_config.generate_tdms_commands():
        tdms_command = test_case[0]
        run_tdms(*tdms_command.generate_tdms_command())

        reference = HDF5File(test_case[1])
        assert HDF5File(tdms_command.output_file).matches(reference)
    return
