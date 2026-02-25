"""
Regression tests for the vectorized smooth-component density loop.

These tests pin the float64 outputs of dmdsm_d2dm from the current scalar
implementation (density_2001_smooth_comps called 500 times per LoS).
After vectorizing with density_2001_smooth_comps_vec the results must match
within a relative tolerance of 1e-7.  A small FP difference is expected
because array vs scalar evaluation of model.arm_radius_splines[j](theta) can
produce results that differ at the ~1e-14 level; these differences propagate
through the CubicSpline interpolation of ne_smooth(s) and ultimately shift DM
and SM by at most ~1e-10, well inside the 1e-7 tolerance.

Cases are identical to test_smallscale_vectorize.py so that any accumulation
of both vectorization changes is covered in either test file.

Cases:
  - in-plane long LoS hitting clumps/voids  (gl=30°,  gb=0°,    d=50 kpc)
  - high-latitude long LoS                  (gl=65°,  gb=10°,   d=50 kpc)
  - moderate distance off-plane             (gl=120°, gb=25°,   d=1.5 kpc)
  - LISM-dominated short LoS (J0323+3944)   (gl=152.18°, gb=-14.338°, d=0.95 kpc)
"""
import pytest
import numpy as np
from numpy import deg2rad
from mwprop.nemod.dmdsm import dmdsm_d2dm

REL = 1e-7    # allow minor FP variation from array vs scalar spline evaluation
ABS = 1e-12

CASES = [
    {
        "name": "inplane_50kpc",
        "ldeg": 30.0, "bdeg": 0.0, "d": 50.0,
        "dm":      1446.7678216162012,
        "sm":      8.5444257036687414,
        "smtau":   6.6315500476504612,
        "smtheta": 18.3433973830609,
        "smiso":   269.37267702696971,
    },
    {
        "name": "offplane_50kpc",
        "ldeg": 65.0, "bdeg": 10.0, "d": 50.0,
        "dm":      134.00308115340957,
        "sm":      4.5439124029864395e-4,
        "smtau":   2.2440369029998707e-4,
        "smtheta": 1.1232530653137534e-3,
        "smiso":   6.9248055402793629e-3,
    },
    {
        "name": "moderate_1p5kpc",
        "ldeg": 120.0, "bdeg": 25.0, "d": 1.5,
        "dm":      20.954337826482909,
        "sm":      8.6709283724333442e-5,
        "smtau":   9.5814025238290896e-5,
        "smtheta": 5.2655340157769143e-5,
        "smiso":   8.1582771369937196e-5,
    },
    {
        "name": "lism_0p95kpc",
        "ldeg": 152.18, "bdeg": -14.338, "d": 0.95,
        "dm":      16.240964985748889,
        "sm":      1.3802585031930291e-4,
        "smtau":   1.5670789856281705e-4,
        "smtheta": 6.3477819939699418e-5,
        "smiso":   6.5942541166867747e-5,
    },
]


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_d2dm_smoothscale_regression(case):
    """
    dmdsm_d2dm must return the same DM and SM values as the scalar
    smooth-component loop within REL=1e-7.
    """
    limit, d_out, dm, sm, smtau, smtheta, smiso = dmdsm_d2dm(
        deg2rad(case["ldeg"]),
        deg2rad(case["bdeg"]),
        case["d"],
        ds_coarse=0.1,
        ds_fine=0.01,
        Nsmin=20,
        d2dm_only=False,
        do_analysis=False,
        plotting=False,
        verbose=False,
    )
    assert dm      == pytest.approx(case["dm"],      rel=REL, abs=ABS)
    assert sm      == pytest.approx(case["sm"],      rel=REL, abs=ABS)
    assert smtau   == pytest.approx(case["smtau"],   rel=REL, abs=ABS)
    assert smtheta == pytest.approx(case["smtheta"], rel=REL, abs=ABS)
    assert smiso   == pytest.approx(case["smiso"],   rel=REL, abs=ABS)
