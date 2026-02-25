import os

import pytest

from mwprop.nemod import config_nemod


def _read_which_model_file():
    params_dir = os.path.join(os.path.dirname(config_nemod.__file__), "params")
    with open(os.path.join(params_dir, "which_model.inp"), "r") as f:
        return f.read().strip()


def test_set_model_switches_flags_and_data():
    # Switch to NE2001 and verify flags/data reflect gal01.inp
    config_nemod.set_model("ne2001")
    assert config_nemod.eval_NE2001 is True
    assert config_nemod.eval_NE2025 is False
    assert _read_which_model_file() == "NE2001"
    assert config_nemod.Dgal["n1h1"] == pytest.approx(0.033)



def test_set_model_propagates_to_density_components():
    from mwprop.nemod import density_components

    config_nemod.set_model("ne2001")
    n1h1_2001 = density_components.n1h1

    config_nemod.set_model("ne2025")
    assert density_components.n1h1 != n1h1_2001
    assert density_components.n1h1 == pytest.approx(config_nemod.n1h1)

    # Switch to NE2025 and verify flags/data reflect gal25.inp
    config_nemod.set_model("ne2025")
    assert config_nemod.eval_NE2025 is True
    assert config_nemod.eval_NE2001 is False
    assert _read_which_model_file() == "NE2025"
    assert config_nemod.Dgal["n1h1"] == pytest.approx(0.0275)

    # Switch back to NE2001 and verify flags/data reflect gal01.inp
    config_nemod.set_model("ne2001")
    assert config_nemod.eval_NE2001 is True
    assert config_nemod.eval_NE2025 is False
    assert _read_which_model_file() == "NE2001"
    assert config_nemod.Dgal["n1h1"] == pytest.approx(0.033)

    # Restore default model for downstream tests
    config_nemod.set_model("ne2025")
