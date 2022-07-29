import pytest
from utils import run_tdms


def test_help_message_prints():

    output = run_tdms("-h", return_output=True)
    assert "help" in output.lower()


@pytest.mark.parametrize("n_args", (1, 4, 5))
def test_exits_with_incorrect_num_arguments(n_args):
    """Test with the incorrect number of arguments"""

    args = n_args * ["a"]
    returncode = run_tdms(*args, return_returncode=True)
    assert returncode != 0
