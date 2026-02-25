#!/usr/bin/env python
"""
Benchmark dmdsm_dm2d with multiple calls to measure amortized JIT speedup.
First run is slow (JIT compilation), subsequent runs are fast (compiled code).
"""

import numpy as np
import time
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from mwprop.nemod.dmdsm import dmdsm_dm2d

# Test cases
test_cases = [
    (np.deg2rad(65.0), np.deg2rad(10.0), 30.0),
    (np.deg2rad(45.0), np.deg2rad(-5.0), 50.0),
    (np.deg2rad(120.0), np.deg2rad(25.0), 20.0),
    (np.deg2rad(200.0), np.deg2rad(-15.0), 40.0),
    (np.deg2rad(350.0), np.deg2rad(8.0), 25.0),
]

print("=" * 80)
print("BENCHMARK: dmdsm_dm2d with Numba JIT")
print("=" * 80)
print()

# First call (includes JIT compilation time)
print("FIRST CALL (includes JIT compilation)...")
t0 = time.time()
l, b, dm = test_cases[0]
limit, dist, dm_calc, sm, smtau, smtheta, smiso = dmdsm_dm2d(
    l, b, dm,
    dm2d_only=False,
    do_analysis=False,
    plotting=False,
    verbose=False,
    debug=False
)
t1 = time.time()
first_call_time = t1 - t0
print(f"  Time: {first_call_time:.4f}s (includes JIT compilation)")
print(f"  Result: d={dist:.2f} kpc, DM={dm_calc:.2f} pc/cm³")
print()

# Subsequent calls (only execution time, JIT already compiled)
print("SUBSEQUENT CALLS (JIT already compiled)...")
times = []
for i, (l, b, dm) in enumerate(test_cases[1:], 1):
    t0 = time.time()
    limit, dist, dm_calc, sm, smtau, smtheta, smiso = dmdsm_dm2d(
        l, b, dm,
        dm2d_only=False,
        do_analysis=False,
        plotting=False,
        verbose=False,
        debug=False
    )
    t1 = time.time()
    elapsed = t1 - t0
    times.append(elapsed)
    print(f"  Call {i}: {elapsed:.4f}s - d={dist:.2f} kpc, DM={dm_calc:.2f}")

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
avg_compiled_time = np.mean(times)
print(f"First call (with JIT): {first_call_time:.4f}s")
print(f"Average compiled call: {avg_compiled_time:.4f}s (n={len(times)})")
print(f"First call overhead:   {first_call_time - avg_compiled_time:.4f}s")
print(f"Speedup ratio:         {first_call_time / avg_compiled_time:.2f}×")
