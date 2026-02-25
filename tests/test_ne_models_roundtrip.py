import pytest

from mwprop.nemod.NE2001 import ne2001
from mwprop.nemod.NE2025 import ne2025


@pytest.mark.parametrize(
    "model_func, expected_dist, expected_dm_roundtrip",
    [
        (ne2001, 3.2587940591703823, 99.99330830834366),
        (ne2025, 2.6198777416235903, 100.00400417623374),
    ],
)
def test_ne_models_roundtrip_dm_distance(model_func, expected_dist, expected_dm_roundtrip):
    ldeg = 200.0
    bdeg = -6.5

    # ndir=1: DM -> D
    _, Dv, _, _ = model_func(
        ldeg,
        bdeg,
        100.0,
        1,
        dmd_only=True,
        classic=False,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=False,
    )

    dist = float(Dv["DIST"])
    dm_out = float(Dv["DM"])

    assert dm_out == pytest.approx(100.0, rel=1e-12, abs=1e-12)
    assert dist == pytest.approx(expected_dist, rel=1e-12, abs=1e-12)

    # ndir=-1: D -> DM using output distance
    _, Dv2, _, _ = model_func(
        ldeg,
        bdeg,
        dist,
        -1,
        dmd_only=True,
        classic=False,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=False,
    )

    dist_roundtrip = float(Dv2["DIST"])
    dm_roundtrip = float(Dv2["DM"])

    assert dist_roundtrip == pytest.approx(dist, rel=1e-12, abs=1e-12)
    assert dm_roundtrip == pytest.approx(expected_dm_roundtrip, rel=1e-12, abs=1e-12)