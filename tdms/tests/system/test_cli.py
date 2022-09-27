import pytest
from utils import run_tdms


def test_help_message_prints():
    """Test that we get the help message"""
    result = run_tdms("-h")
    assert "help" in result.stdout.lower()


@pytest.mark.parametrize("n_args", (1, 4, 5))
def test_exits_with_incorrect_num_arguments(n_args):
    """Test with the incorrect number of arguments"""
    args = n_args * ["a"]
    result = run_tdms(*args)
    assert result.return_code != 0
