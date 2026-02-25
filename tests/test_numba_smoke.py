"""Smoke tests for numba njit functions (run w/o numba to check coverage)"""
import pytest


def test_numba_smoke_nevoidN():
    pytest.importorskip("numba")

    from mwprop.nemod import config_nemod
    from mwprop.nemod.nevoidN import nevoidN, _nevoidN_jit

    if config_nemod.nvoids == 0:
        pytest.skip("No voids defined in model parameters.")

    x = config_nemod.xv[0]
    y = config_nemod.yv[0]
    z = config_nemod.zv[0]

    nevoidN(x, y, z)
    assert _nevoidN_jit.signatures, "Expected numba to compile _nevoidN_jit"


def test_numba_smoke_density_components():
    pytest.importorskip("numba")

    from mwprop.nemod.density_components import ne_outer, ne_inner, _ne_outer_jit, _ne_inner_jit

    ne_outer(0.1, 0.2, 0.0)
    ne_inner(0.1, 0.2, 0.0)

    assert _ne_outer_jit.signatures, "Expected numba to compile _ne_outer_jit"
    assert _ne_inner_jit.signatures, "Expected numba to compile _ne_inner_jit"