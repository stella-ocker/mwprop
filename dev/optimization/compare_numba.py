#!/usr/bin/env python
"""
Compare: Numba JIT vs Pure Python
This tests if Numba actually speeds up the execution (once compiled).
"""

import numpy as np
import time
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

print("=" * 80)
print("COMPARISON: Numba JIT vs Pure Python (compiled calls only)")
print("=" * 80)
print()

# Test cases  
test_cases = [
    (np.deg2rad(65.0), np.deg2rad(10.0), 30.0),
    (np.deg2rad(45.0), np.deg2rad(-5.0), 50.0),
    (np.deg2rad(120.0), np.deg2rad(25.0), 20.0),
    (np.deg2rad(200.0), np.deg2rad(-15.0), 40.0),
    (np.deg2rad(350.0), np.deg2rad(8.0), 25.0),
    (np.deg2rad(30.0), np.deg2rad(12.0), 35.0),
    (np.deg2rad(150.0), np.deg2rad(-20.0), 45.0),
    (np.deg2rad(270.0), np.deg2rad(5.0), 22.0),
]

from mwprop.nemod.dmdsm import dmdsm_dm2d

print("Warming up JIT compilation (first calls)...")
for i, (l, b, dm) in enumerate(test_cases[:2]):
    dmdsm_dm2d(l, b, dm, dm2d_only=False, do_analysis=False, 
               plotting=False, verbose=False, debug=False)
print("JIT compiled.\n")

print("RUNNING BENCHMARK (6 iterations, 8 test cases each)...")
times_numba = []
for iteration in range(6):
    iter_times = []
    for l, b, dm in test_cases:
        t0 = time.perf_counter()
        limit, dist, dm_calc, sm, smtau, smtheta, smiso = dmdsm_dm2d(
            l, b, dm, 
            dm2d_only=False, 
            do_analysis=False, 
            plotting=False, 
            verbose=False, 
            debug=False
        )
        t1 = time.perf_counter()
        iter_times.append(t1 - t0)
    times_numba.append(iter_times)

times_numba = np.array(times_numba)

print()
print("=" * 80)
print("RESULTS (Numba JIT - Compiled Calls Only)")
print("=" * 80)
print(f"Total time over {len(test_cases)} cases × 6 iterations: {times_numba.sum():.4f}s")
print(f"Average per call: {times_numba.mean():.4f}s")
print(f"Min per call: {times_numba.min():.4f}s")
print(f"Max per call: {times_numba.max():.4f}s")
print(f"Std dev: {times_numba.std():.4f}s")
print()
print("Iteration breakdown:")
for i, iter_times in enumerate(times_numba):
    print(f"  Iteration {i+1}: avg={iter_times.mean():.4f}s, "
          f"min={iter_times.min():.4f}s, max={iter_times.max():.4f}s")
