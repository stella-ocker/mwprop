"""
Baseline image tests for matplotlib plots.
"""
import pytest
from numpy import deg2rad

from mwprop.nemod.dmdsm import dmdsm_dm2d, plot_dm_along_LoS
from mwprop.nemod.neclumpN_fast import relevant_clumps
from mwprop.nemod.config_nemod import rcmult, dc


@pytest.mark.mpl_image_compare(tolerance=20)
def test_dm2d_plot_baseline():
    """
    Baseline image for DM vs distance plot.
    """
    l = deg2rad(30.0)
    b = deg2rad(0.0)
    dm_target = 100.0

    limit, dhat, dm_return, sf_vec, dm_cumulate_vec = dmdsm_dm2d(
        l,
        b,
        dm_target,
        ds_fine=0.05,
        ds_coarse=0.2,
        Nsmin=20,
        dm2d_only=True,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=True,
    )

    relevant_clump_indices = relevant_clumps(l, b, sf_vec[-1], rcmult)

    fig = plot_dm_along_LoS(
        dm_target,
        dhat,
        sf_vec,
        dm_cumulate_vec,
        relevant_clump_indices,
        dc,
        which="dm2d",
        plot_dm_target=True,
        saveplot=False,
        show_plot=False,
        annotate_stamp=False,
        return_fig=True,
    )

    return fig
