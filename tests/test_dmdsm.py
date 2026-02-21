"""
Tests for dmdsm.py module
Based on the original __main__ code
"""
import pytest
import numpy as np
from numpy import deg2rad, rad2deg
from mwprop.nemod.dmdsm import dmdsm_dm2d


@pytest.mark.plot
def test_dmdsm_dm2d_basic():
    """
    Test dmdsm_dm2d function with basic parameters
    Original test case from __main__: ldeg=30, bdeg=0, dm_target=1000
    """
    # Test parameters from original __main__
    ldeg, bdeg, dm_target = 30, 0, 1000
    ndir = -1
    ds_fine = 0.005
    ds_coarse = 0.1

    l = deg2rad(ldeg)
    b = deg2rad(bdeg)

    verbose = True
    do_analysis = True
    dm2d_only = False
    plotting = True
    debug = False

    # Run the function
    limit, dhat, dm_target_out, sm, smtau, smtheta, smiso = dmdsm_dm2d(
        l, b, dm_target,
        ds_fine=ds_fine,
        ds_coarse=ds_coarse,
        Nsmin=10,
        dm2d_only=dm2d_only,
        do_analysis=do_analysis,
        plotting=plotting,
        verbose=verbose,
        debug=debug
    )

    # Basic assertions to verify output
    # Note: limit can be a single character or string (possibly with whitespace)
    assert limit is not None, "limit should not be None"
    assert isinstance(limit, str), "limit should be a string"
    assert dhat > 0, f"Distance should be positive, got {dhat}"
    assert dm_target_out == dm_target, f"DM target mismatch: {dm_target_out} != {dm_target}"
    assert sm >= 0, f"SM should be non-negative, got {sm}"
    assert smtau >= 0, f"SMtau should be non-negative, got {smtau}"
    assert smtheta >= 0, f"SMtheta should be non-negative, got {smtheta}"
    assert smiso >= 0, f"SMiso should be non-negative, got {smiso}"

    print(f"\nTest results:")
    print(f"  limit: {limit}")
    print(f"  dhat: {dhat:.3f} kpc")
    print(f"  DM: {dm_target_out:.1f} pc/cc")
    print(f"  SM: {sm:.4f}")
    print(f"  SMtau: {smtau:.4f}")
    print(f"  SMtheta: {smtheta:.4f}")
    print(f"  SMiso: {smiso:.4f}")


@pytest.mark.plot
def test_dmdsm_dm2d_only_mode():
    """
    Test dmdsm_dm2d function in dm2d_only mode (faster, no scattering)
    """
    ldeg, bdeg, dm_target = 30, 0, 1000
    ds_fine = 0.005
    ds_coarse = 0.1

    l = deg2rad(ldeg)
    b = deg2rad(bdeg)

    verbose = True
    do_analysis = True
    dm2d_only = True
    plotting = True
    debug = False

    # Run the function in dm2d_only mode
    limit, dhat, dm_target_out = dmdsm_dm2d(
        l, b, dm_target,
        ds_fine=ds_fine,
        ds_coarse=ds_coarse,
        Nsmin=10,
        dm2d_only=dm2d_only,
        do_analysis=do_analysis,
        plotting=plotting,
        verbose=verbose,
        debug=debug
    )

    # Basic assertions
    assert limit is not None, "limit should not be None"
    assert isinstance(limit, str), "limit should be a string"
    assert dhat > 0, f"Distance should be positive, got {dhat}"
    assert dm_target_out == dm_target, f"DM target mismatch: {dm_target_out} != {dm_target}"

    print(f"\nTest results (dm2d_only mode):")
    print(f"  limit: {limit}")
    print(f"  dhat: {dhat:.3f} kpc")
    print(f"  DM: {dm_target_out:.1f} pc/cc")
