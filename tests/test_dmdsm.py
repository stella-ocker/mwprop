"""
Tests for dmdsm.py module
Based on the original __main__ code
"""
import pytest
import numpy as np
from numpy import deg2rad, rad2deg
from mwprop.nemod.dmdsm import dmdsm_dm2d


DM2D_ONLY_DATASETS = [
    {
        "name": "case_30_0_100",
        "ldeg": 30.0,
        "bdeg": 0.0,
        "dm": 100.0,
        "dhat": 2.948350355206938,
        "limit": " ",
    },
    {
        "name": "case_45_5_200",
        "ldeg": 45.0,
        "bdeg": 5.0,
        "dm": 200.0,
        "dhat": 12.758739084001773,
        "limit": " ",
    },
    {
        "name": "case_120_-2_50",
        "ldeg": 120.0,
        "bdeg": -2.0,
        "dm": 50.0,
        "dhat": 2.541696720722774,
        "limit": " ",
    },
]

FULL_OUTPUT_DATASET = {
    "name": "case_30_0_100_full",
    "ldeg": 30.0,
    "bdeg": 0.0,
    "dm": 100.0,
    "dhat": 2.948350355206938,
    "limit": " ",
    "sm": 0.09242970182344527,
    "smtau": 0.06086291443591926,
    "smtheta": 0.01028680389779793,
    "smiso": 0.43524360221382075,
}


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


@pytest.mark.parametrize(
    "dataset",
    DM2D_ONLY_DATASETS,
    ids=[d["name"] for d in DM2D_ONLY_DATASETS],
)
def test_dmdsm_dm2d_only_expected(dataset):
    """
    Regression tests for dm2d_only outputs against precomputed datasets.
    """
    l = deg2rad(dataset["ldeg"])
    b = deg2rad(dataset["bdeg"])

    limit, dhat, dm_target_out = dmdsm_dm2d(
        l,
        b,
        dataset["dm"],
        ds_fine=0.05,
        ds_coarse=0.2,
        Nsmin=20,
        dm2d_only=True,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=False,
    )

    assert limit == dataset["limit"]
    assert dm_target_out == dataset["dm"]
    assert dhat == pytest.approx(dataset["dhat"], rel=1e-4, abs=1e-6)


def test_dmdsm_full_outputs_expected():
    """
    Regression test for full dmdsm outputs (DM + SM variants).
    """
    l = deg2rad(FULL_OUTPUT_DATASET["ldeg"])
    b = deg2rad(FULL_OUTPUT_DATASET["bdeg"])

    limit, dhat, dm_target_out, sm, smtau, smtheta, smiso = dmdsm_dm2d(
        l,
        b,
        FULL_OUTPUT_DATASET["dm"],
        ds_fine=0.05,
        ds_coarse=0.2,
        Nsmin=20,
        dm2d_only=False,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=False,
    )

    assert limit == FULL_OUTPUT_DATASET["limit"]
    assert dm_target_out == FULL_OUTPUT_DATASET["dm"]
    assert dhat == pytest.approx(FULL_OUTPUT_DATASET["dhat"], rel=1e-4, abs=1e-6)
    assert sm == pytest.approx(FULL_OUTPUT_DATASET["sm"], rel=1e-4, abs=1e-6)
    assert smtau == pytest.approx(FULL_OUTPUT_DATASET["smtau"], rel=1e-4, abs=1e-6)
    assert smtheta == pytest.approx(FULL_OUTPUT_DATASET["smtheta"], rel=1e-4, abs=1e-6)
    assert smiso == pytest.approx(FULL_OUTPUT_DATASET["smiso"], rel=1e-4, abs=1e-6)
