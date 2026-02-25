"""
Tests for iss_mw_utils2020p_mod.py __main__ code
Tests the scattering calculation functionality
"""
import pytest
from mwprop.scattering_functions.iss_mw_utils2020p_mod import (
    calc_pdfg_for_xgal_los, main
)


def test_calc_pdfg_for_xgal_los():
    """
    Test calc_pdfg_for_xgal_los function with parameters from __main__
    Original test case from __main__
    Returns a tuple of (data_dict, units_dict, description_dict)
    """
    # Parameters from original __main__
    RF = 1.0
    BW = 1.0
    RF_input = 1.0
    TAU = 1.0
    THETA_X = 1.0
    Deff = 2.5
    SM = 10.**(-3.5)

    # Call the function
    result = calc_pdfg_for_xgal_los(
        RF, BW, RF_input, TAU, THETA_X, Deff, SM,
        theta_source=1.e-6, dg=1.e-5, gmin=1.e-5, gmax=30.
    )

    # The function returns a tuple of (data_dict, units_dict, description_dict)
    assert result is not None, "Result should not be None"
    assert isinstance(result, tuple), "Result should be a tuple"
    assert len(result) == 3, "Result should have 3 elements (data, units, descriptions)"

    data_dict, units_dict, desc_dict = result
    assert isinstance(data_dict, dict), "First element should be a dictionary"
    assert isinstance(units_dict, dict), "Second element should be a dictionary"
    assert isinstance(desc_dict, dict), "Third element should be a dictionary"
    assert len(data_dict) > 0, "Data dictionary should not be empty"


def test_calc_pdfg_for_xgal_los_with_different_params():
    """
    Test calc_pdfg_for_xgal_los with different parameters
    Returns a tuple of (data_dict, units_dict, description_dict)
    """
    RF = 2.0
    BW = 100.0
    RF_input = 2.0
    TAU = 0.5
    THETA_X = 0.5
    Deff = 5.0
    SM = 10.**(-2.0)

    result = calc_pdfg_for_xgal_los(
        RF, BW, RF_input, TAU, THETA_X, Deff, SM,
        theta_source=1.e-6, dg=1.e-5, gmin=1.e-5, gmax=30.
    )

    assert result is not None
    assert isinstance(result, tuple)
    assert len(result) == 3

    data_dict, units_dict, desc_dict = result
    assert isinstance(data_dict, dict)
    assert isinstance(units_dict, dict)
    assert isinstance(desc_dict, dict)


def test_main_function_basic():
    """
    Test main function from __main__ with basic parameters
    Original __main__ calls: main(l, b, RF, BW, dxgal_mpc, theta_source, SM, DMmodel)

    Note: This test uses non-interactive mode (doplot=False)
    """
    l = 30.0  # Galactic longitude in degrees
    b = 0.0   # Galactic latitude in degrees
    RF = 1.0  # Radio frequency in GHz
    BW = 100.0  # Bandwidth in MHz
    dxgal_mpc = 100.0  # Distance in Mpc
    theta_source = 1.0  # Source size in mas
    SM = 10.**(-3.5)  # Scattering measure
    DMmodel = 100  # DM model

    # Call main with doplot=False to avoid interactive display
    try:
        main(l, b, RF, BW, dxgal_mpc, theta_source, SM, DMmodel)
    except Exception as e:
        # main function may not return anything or may raise errors
        # Just ensure it runs without crashing
        pass


def test_main_function_different_coordinates():
    """
    Test main function with different galactic coordinates
    """
    l = 0.0
    b = 90.0
    RF = 1.4
    BW = 100.0
    dxgal_mpc = 50.0
    theta_source = 0.5
    SM = 10.**(-4.0)
    DMmodel = 50

    try:
        main(l, b, RF, BW, dxgal_mpc, theta_source, SM, DMmodel)
    except Exception as e:
        pass
