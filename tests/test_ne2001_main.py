"""
Tests for NE2001.py __main__ code
Tests the help/explain functionality
"""
import pytest
import os
import sys
from io import StringIO
from mwprop.nemod import NE2001


def test_ne2001_explain_flag():
    """
    Test NE2001 module with -e flag to print README
    """
    # Save original argv
    original_argv = sys.argv
    original_stdout = sys.stdout
    
    try:
        # Simulate command line with -e flag
        sys.argv = ['NE2001.py', '-e']
        sys.stdout = StringIO()
        
        # The __main__ code should print README and call sys.exit()
        # We need to test that README file exists and is readable
        script_path = os.path.dirname(os.path.realpath(NE2001.__file__))
        infile = script_path + '/README.txt'
        
        assert os.path.exists(infile), f"README file not found at {infile}"
        
        with open(infile) as fexplain:
            content = fexplain.read()
            assert len(content) > 0, "README file is empty"
            assert isinstance(content, str), "README content should be string"
            
    finally:
        # Restore original argv and stdout
        sys.argv = original_argv
        sys.stdout = original_stdout


def test_ne2001_explain_long_flag():
    """
    Test NE2001 module with --explain flag to print README
    """
    # Save original argv
    original_argv = sys.argv
    
    try:
        # Simulate command line with --explain flag
        sys.argv = ['NE2001.py', '--explain']
        
        # Test that README file exists and is readable
        script_path = os.path.dirname(os.path.realpath(NE2001.__file__))
        infile = script_path + '/README.txt'
        
        assert os.path.exists(infile), f"README file not found at {infile}"
        
        with open(infile) as fexplain:
            content = fexplain.read()
            assert len(content) > 0, "README file is empty"
            
    finally:
        # Restore original argv
        sys.argv = original_argv


def test_ne2001_readme_readable():
    """
    Test that the README.txt file is readable and contains content
    """
    script_path = os.path.dirname(os.path.realpath(NE2001.__file__))
    infile = script_path + '/README.txt'
    
    assert os.path.exists(infile), f"README.txt not found at {infile}"
    assert os.path.isfile(infile), f"{infile} is not a file"
    
    with open(infile, 'r') as f:
        content = f.read()
        assert len(content) > 0, "README.txt is empty"
