"""
Tests for ne_input.py __main__ code
Tests the read_nemod_parameters function
"""
import pytest
from mwprop.nemod.ne_input import read_nemod_parameters


def test_read_nemod_parameters():
    """
    Test read_nemod_parameters function
    Original __main__ code: Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()
    """
    # Call the function as in the original __main__
    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()
    
    # Verify that we got 8 return values
    assert Dgal is not None, "Dgal should not be None"
    assert Dgc is not None, "Dgc should not be None"
    assert Dlism is not None, "Dlism should not be None"
    assert Dclumps is not None, "Dclumps should not be None"
    assert Dvoids is not None, "Dvoids should not be None"
    assert Darms is not None, "Darms should not be None"
    assert eval_NE2025 is not None, "eval_NE2025 should not be None"
    assert eval_NE2001 is not None, "eval_NE2001 should not be None"


def test_read_nemod_parameters_returns_dicts():
    """
    Test that read_nemod_parameters returns dictionary-like objects
    """
    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()
    
    # Check that the first 6 returns are dictionaries
    assert isinstance(Dgal, dict), "Dgal should be a dictionary"
    assert isinstance(Dgc, dict), "Dgc should be a dictionary"
    assert isinstance(Dlism, dict), "Dlism should be a dictionary"
    assert isinstance(Dclumps, dict), "Dclumps should be a dictionary"
    assert isinstance(Dvoids, dict), "Dvoids should be a dictionary"
    assert isinstance(Darms, dict), "Darms should be a dictionary"


def test_read_nemod_parameters_eval_flags():
    """
    Test that eval_NE2025 and eval_NE2001 are boolean flags
    """
    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()
    
    # Check that the last 2 returns are booleans
    assert isinstance(eval_NE2025, (bool, int)), "eval_NE2025 should be a boolean or int"
    assert isinstance(eval_NE2001, (bool, int)), "eval_NE2001 should be a boolean or int"


def test_read_nemod_parameters_dicts_have_content():
    """
    Test that returned dictionaries contain keys
    """
    Dgal, Dgc, Dlism, Dclumps, Dvoids, Darms, eval_NE2025, eval_NE2001 = read_nemod_parameters()
    
    # Check that dictionaries are not empty
    assert len(Dgal) > 0, "Dgal dictionary should not be empty"
    assert len(Dgc) > 0, "Dgc dictionary should not be empty"
    assert len(Dlism) > 0, "Dlism dictionary should not be empty"
    assert len(Dclumps) > 0, "Dclumps dictionary should not be empty"
    assert len(Dvoids) > 0, "Dvoids dictionary should not be empty"
    assert len(Darms) > 0, "Darms dictionary should not be empty"
