import pytest

from mwprop.nemod import config_nemod
from mwprop.nemod.nevoidN import nevoidN


def test_nevoidN_hits_void_center():
    if config_nemod.nvoids == 0:
        pytest.skip("No voids defined in model parameters.")

    x = config_nemod.xv[0]
    y = config_nemod.yv[0]
    z = config_nemod.zv[0]

    nevN, FvN, hitvoid, wvoid = nevoidN(x, y, z)

    assert hitvoid == 1
    assert wvoid == 1
    assert nevN == pytest.approx(config_nemod.nev[0])
    assert FvN == pytest.approx(config_nemod.Fv[0])