import pytest

from mwprop.nemod import config_nemod
from mwprop.nemod.density_components import ne_gc


def test_ne_gc_at_center():
    ne_gc_out, F_gc = ne_gc(config_nemod.xgc, config_nemod.ygc, config_nemod.zgc)

    assert ne_gc_out == pytest.approx(config_nemod.negc0)
    assert F_gc == pytest.approx(config_nemod.Fgc0)


def test_ne_gc_far_from_center():
    far_y = config_nemod.ygc + 5.0 * config_nemod.rgc
    ne_gc_out, F_gc = ne_gc(config_nemod.xgc, far_y, config_nemod.zgc)

    assert ne_gc_out == 0
    assert F_gc == 0