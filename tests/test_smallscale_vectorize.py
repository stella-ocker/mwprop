"""
Regression tests for the vectorized smallscale density loop.

These tests pin the exact float64 outputs of dmdsm_d2dm BEFORE the
vectorization refactor.  After the refactor the results must be bitwise
identical (same floating-point operations, same order), so we use a very
tight relative tolerance of 1e-10.

Cases are chosen to exercise:
  - in-plane long LoS that hits clumps and voids  (gl=30, gb=0, d=50 kpc)
  - high-latitude long LoS                         (gl=65, gb=10, d=50 kpc)
  - moderate distance off-plane                    (gl=120, gb=25, d=1.5 kpc)
  - LISM-dominated short LoS (J0323+3944 direction) (gl=152.18, gb=-14.338, d=0.95 kpc)
"""
import pytest
import numpy as np
from numpy import deg2rad
from mwprop.nemod.dmdsm import dmdsm_d2dm

REL = 1e-10   # tight: same arithmetic must give same bits
ABS = 1e-15

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
def test_d2dm_smallscale_regression(case):
    """
    dmdsm_d2dm must return exactly the same DM and SM values as before
    the smallscale-loop vectorization.
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

    assert limit == " "
    assert float(dm)     == pytest.approx(case["dm"],      rel=REL, abs=ABS)
    assert float(sm)     == pytest.approx(case["sm"],      rel=REL, abs=ABS)
    assert float(smtau)  == pytest.approx(case["smtau"],   rel=REL, abs=ABS)
    assert float(smtheta)== pytest.approx(case["smtheta"], rel=REL, abs=ABS)
    assert float(smiso)  == pytest.approx(case["smiso"],   rel=REL, abs=ABS)
