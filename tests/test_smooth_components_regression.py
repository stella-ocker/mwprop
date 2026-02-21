"""
Regression tests for smooth-component helpers prior to performance optimizations.
"""
import pytest

from mwprop.nemod.config_nemod import Ncoarse
from mwprop.nemod.ne_arms import ne_arms_ne2001p
from mwprop.nemod.density import density_2001_smooth_comps


NE_ARMS_DATASETS = [
    {
        "point": (0.0, 5.0, 0.0),
        "nea": 0.0009229400428876469,
        "F": 3.7,
        "whicharm": 2,
    },
    {
        "point": (2.0, 2.0, 0.1),
        "nea": 0.0011933057859719425,
        "F": 3.7,
        "whicharm": 2,
    },
    {
        "point": (-3.0, 5.0, -0.2),
        "nea": 0.014112196847313123,
        "F": 3.7,
        "whicharm": 2,
    },
    {
        "point": (1.5, 9.0, 0.3),
        "nea": 0.0037533252616640065,
        "F": 3.7,
        "whicharm": 5,
    },
]


SMOOTH_COMP_DATASETS = [
    {
        "point": (0.0, 5.0, 0.0),
        "ne1": 0.02157265129165555,
        "ne2": 0.06877168237140906,
        "nea": 0.0009229400428876469,
        "F1": 0.18,
        "F2": 120.0,
        "Fa": 3.7,
        "whicharm": 2,
        "ne_smooth": 0.09126727370595225,
        "Fsmooth": 68.14545497955861,
    },
    {
        "point": (2.0, 2.0, 0.1),
        "ne1": 0.023084779094999033,
        "ne2": 0.027077149102683155,
        "nea": 0.0011933057859719425,
        "F1": 0.18,
        "F2": 120.0,
        "Fa": 3.7,
        "whicharm": 2,
        "ne_smooth": 0.05135523398365413,
        "Fsmooth": 33.39772745791793,
    },
    {
        "point": (-3.0, 5.0, -0.2),
        "ne1": 0.02041342683380049,
        "ne2": 0.009430069677918667,
        "nea": 0.014112196847313123,
        "F1": 0.18,
        "F2": 120.0,
        "Fa": 3.7,
        "whicharm": 2,
        "ne_smooth": 0.04395569335903228,
        "Fsmooth": 5.943277056652071,
    },
    {
        "point": (1.5, 9.0, 0.3),
        "ne1": 0.015783505205136154,
        "ne2": 4.292838725444624e-06,
        "nea": 0.0037533252616640065,
        "F1": 0.18,
        "F2": 120.0,
        "Fa": 3.7,
        "whicharm": 5,
        "ne_smooth": 0.019541123305525605,
        "Fsmooth": 0.25393690783409784,
    },
]


@pytest.mark.parametrize(
    "dataset",
    NE_ARMS_DATASETS,
    ids=[f"ne_arms_{i}" for i in range(len(NE_ARMS_DATASETS))],
)
def test_ne_arms_ne2001p_regression(dataset):
    x, y, z = dataset["point"]
    nea, F, whicharm = ne_arms_ne2001p(x, y, z, Ncoarse=Ncoarse)

    assert nea == pytest.approx(dataset["nea"], rel=1e-12, abs=1e-12)
    assert F == pytest.approx(dataset["F"], rel=1e-12, abs=1e-12)
    assert int(whicharm) == dataset["whicharm"]


@pytest.mark.parametrize(
    "dataset",
    SMOOTH_COMP_DATASETS,
    ids=[f"smooth_comps_{i}" for i in range(len(SMOOTH_COMP_DATASETS))],
)
def test_density_2001_smooth_comps_regression(dataset):
    x, y, z = dataset["point"]
    ne1, ne2, nea, F1, F2, Fa, whicharm, ne_smooth, Fsmooth = density_2001_smooth_comps(x, y, z)

    assert ne1 == pytest.approx(dataset["ne1"], rel=1e-12, abs=1e-12)
    assert ne2 == pytest.approx(dataset["ne2"], rel=1e-12, abs=1e-12)
    assert nea == pytest.approx(dataset["nea"], rel=1e-12, abs=1e-12)
    assert F1 == pytest.approx(dataset["F1"], rel=1e-12, abs=1e-12)
    assert F2 == pytest.approx(dataset["F2"], rel=1e-12, abs=1e-12)
    assert Fa == pytest.approx(dataset["Fa"], rel=1e-12, abs=1e-12)
    assert int(whicharm) == dataset["whicharm"]
    assert ne_smooth == pytest.approx(dataset["ne_smooth"], rel=1e-12, abs=1e-12)
    assert Fsmooth == pytest.approx(dataset["Fsmooth"], rel=1e-12, abs=1e-12)
