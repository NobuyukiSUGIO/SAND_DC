#!/usr/bin/env python3
"""
17-round SAND-64 key recovery: threshold, retained classes and success
probability of the counting phase (Section 7.3.2 and Table 12 of the paper).

Model (Section 7.3.1): N = 2^62 pairs, p = 2^-57, b = 64, k = 69,

    c_*     ~ Poisson(2^5)
    c_bg    ~ Poisson(2^-2)                for 2^69 - 1 background classes
    c_near  ~ Poisson(2^5 rho + 2^-2)      for the nearby classes,

the nearby term being read from nearby/results/nearby_sand64_seed*.json.
Candidate master keys are verified against two known plaintext-ciphertext
pairs: with one 64-bit pair, (2^59 - 1) 2^-64 ~ 2^-5 wrong completions of each
retained class would pass.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.counting import load_nearby, poisson_ge, print_table, threshold_table

LOG_N, LOG_P, BLOCK, K = 62, -57, 64, 69
TAU = 16
TAUS = [10, 12, 14, 15, 16]
ROUNDS = 17
# counting phase: Steps 1, 2 and the three stages (Section 7.3.2)
T_COUNT = math.log2(2 / ROUNDS * (2.0 ** 63 + 2.0 ** 86 + 2.0 ** 85 + 2.0 ** 76))

nearby = load_nearby(64)
if not nearby:
    sys.exit("no nearby/results/nearby_sand64_seed*.json; run nearby/nearby_classes.py 64")
print(f"nearby-class runs: seeds {sorted(nearby)}")
for seed, v in sorted(nearby.items()):
    print(f"  seed {seed}: E[nearby retained] at tau={TAU} is {v[TAU]:.1f}")
ref_seed = 23 if 23 in nearby else min(nearby)

print(f"\nTable 12 (seed {ref_seed}); counting floor = 2^{T_COUNT:.2f}")
rows = threshold_table(LOG_N, LOG_P, BLOCK, K, TAUS, nearby[ref_seed], T_COUNT)
print_table(rows, "T")

bg = {t: (2.0 ** K - 1) * poisson_ge(2.0 ** (LOG_N - BLOCK), t) for t in TAUS}
print(f"\nbackground alone drops below one at tau = {min(t for t in TAUS if bg[t] < 1)}"
      f" ({bg[15]:.3f}); at tau = 16 it is {bg[16]:.4f}")
succ = poisson_ge(2.0 ** (LOG_N + LOG_P), TAU)
n_ret = 1 + bg[TAU] + nearby[ref_seed][TAU]
print(f"tau = {TAU}: Pr[c_* >= tau] = {succ:.5f}; about {n_ret:.0f} classes searched, "
      f"final stage 2^{math.log2(n_ret * 2.0 ** (128 - K)):.1f}")
print(f"verification: one pair -> {(2**59 - 1) * 2.0**-64:.4f} false keys per class, "
      f"two pairs -> 2^{59 - 128}")
print(f"ignoring the clustering (p = 2^-58): Pr[c_* >= {TAU}] = "
      f"{poisson_ge(2.0 ** (LOG_N - 58), TAU):.3f}")
