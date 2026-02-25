#!/usr/bin/env python
"""
Profile dmdsm_dm2d to identify bottlenecks.
Measures time spent in each major density function.
"""

import numpy as np
import cProfile
import pstats
from io import StringIO
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from mwprop.nemod.dmdsm import dmdsm_dm2d
from mwprop.nemod.density import density_2001_smooth_comps, density_2001_smallscale_comps
from mwprop.nemod.ne_arms import ne_arms_ne2001p
from mwprop.nemod.density_components import ne_outer, ne_inner

# Test case
l = np.deg2rad(65.0)   # Galactic longitude
b = np.deg2rad(10.0)   # Galactic latitude
dm_target = 30.0       # pc/cm^3

print("=" * 80)
print("PROFILING dmdsm_dm2d")
print("=" * 80)
print(f"Input: l={np.rad2deg(l):.1f}°, b={np.rad2deg(b):.1f}°, DM={dm_target:.1f} pc/cm³")
print()

# Profile the full call
pr = cProfile.Profile()
pr.enable()

limit, dist, dm_calc, sm, smtau, smtheta, smiso = dmdsm_dm2d(
    l, b, dm_target,
    dm2d_only=False,
    do_analysis=False,
    plotting=False,
    verbose=False,
    debug=False
)

pr.disable()

# Print stats
s = StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
ps.print_stats(30)  # Top 30 functions
print(s.getvalue())

print("\n" + "=" * 80)
print("RESULT")
print("=" * 80)
print(f"Distance: {dist:.2f} kpc")
print(f"DM: {dm_calc:.2f} pc/cm³")
print(f"SM: {sm:.6e}")
