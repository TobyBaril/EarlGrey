"""
To run the tests, use the commands:
export PYTHONPATH=$PYTHONPATH:/Users/feurtey/Local_documents/Projects/EarlGrey_pangenome/
pytest pangenome/tests/on_start_functions_tests.py
"""

import pytest
from pangenome.scripts.on_start_functions import (
    convert_seconds,
    validate_parameters
)


# Test cases for the convert_seconds functions
def test_convert_seconds_zero():
    assert convert_seconds(0) == "00:00:00.00"

def test_convert_seconds_under_minute():
    assert convert_seconds(42) == "00:00:42.00"

def test_convert_seconds_minutes_and_seconds():
    assert convert_seconds(65) == "00:01:05.00"

def test_convert_seconds_hours():
    assert convert_seconds(3661) == "01:01:01.00"

def test_convert_seconds_float():
    assert convert_seconds(12.9) == "00:00:12.00"

def test_convert_seconds_numeric_string():
    assert convert_seconds("3600") == "01:00:00.00"

def test_convert_seconds_invalid_input():
    with pytest.raises(ValueError):
        convert_seconds("tea")


# Test cases for the validate_parameters function
# ---------
def minimal_valid_config():
    return {
        "genome": "genome.fa",
        "species": "human",
        "directory": "/tmp"
    }

# Required parameters tests
def test_validate_parameters_missing_required_exits(capsys):
    """
    Missing required parameter leads to exit
    
    :param capsys: Description
    """
    config = {
        "genome": "genome.fa",
        "species": "human"
        # directory missing
    }

    with pytest.raises(SystemExit) as exc:
        validate_parameters(config)

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "Required parameter 'directory'" in captured.out

def test_validate_parameters_required_present():
    """ 
    All required parameters present passes validation
    """
    config = minimal_valid_config()
    result = validate_parameters(config)
    assert result is config

# Optional parameters tests
def test_validate_parameters_sets_defaults():
    """
    If not provided, optional parameters are set to default values
    If provided, they are preserved
    """
    # Not provided
    config = minimal_valid_config()

    result = validate_parameters(config)

    assert result["ProcNum"] == 1
    assert result["num"] == 10
    assert result["Flank"] == 1000
    assert result["min_seq"] == 3

    # Provided
    config = minimal_valid_config()
    config["ProcNum"] = 8
    config["cluster"] = "yes"

    result = validate_parameters(config)

    assert result["ProcNum"] == 8
    assert result["cluster"] == "yes"


def test_validate_parameters_repeatmasker(capsys):
    """
    Check that the appropriate message is printed based on RepSpec presence
    """
    #No RepSpec
    config = minimal_valid_config()

    validate_parameters(config)
    captured = capsys.readouterr()

    assert "without an initial mask" in captured.out

    #With RepSpec
    config = minimal_valid_config()
    config["RepSpec"] = "human"

    validate_parameters(config)
    captured = capsys.readouterr()

    assert "with an initial mask" in captured.out


def test_validate_parameters_cluster_yes(capsys):
    config = minimal_valid_config()
    config["cluster"] = "yes"

    validate_parameters(config)
    captured = capsys.readouterr()

    assert "will be clustered" in captured.out

def test_validate_parameters_heli_default(capsys):
    config = minimal_valid_config()

    validate_parameters(config)
    captured = capsys.readouterr()

    assert "HELITRON detection will not be run" in captured.out
