#!/usr/bin/env python3
"""
23-round SAND-128 key recovery: threshold, retained classes and success
probability of the counting phase (Section 7.3.3 and Table 13 of the paper).

Model (Section 7.3.1): N = 2^125 pairs, p = 2^-120, b = 128, k = 72,

    c_*     ~ Poisson(2^5)
    c_bg    ~ Poisson(2^-3)                for 2^72 - 1 background classes
    c_near  ~ Poisson(2^5 rho + 2^-3)      for the nearby classes,

the nearby term being read from nearby/results/nearby_sand128_seed*.json
(produced by nearby/nearby_classes.py).  Every retained class is searched, so
the final stage costs (1 + E[N_wrong]) 2^56 and the success probability is
Pr[c_* >= tau] minus a verification false positive of about 2^-72 per class.

This script only re-evaluates formulas; the nearby term is the one measured
quantity it uses.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.counting import load_nearby, poisson_ge, print_table, threshold_table

LOG_N, LOG_P, BLOCK, K = 125, -120, 128, 72
TAU = 14
TAUS = [10, 12, 14, 16]
ROUNDS = 23
# the two stages of Step 4 (Section 7.3.3): 2^35 x 2^46 and 2^72 x 2^9 round pairs
T_STAGES = math.log2(2 / ROUNDS * (2.0 ** 81 + 2.0 ** 81))

nearby = load_nearby(128)
if not nearby:
    sys.exit("no nearby/results/nearby_sand128_seed*.json; run nearby/nearby_classes.py 128")
print(f"nearby-class runs: seeds {sorted(nearby)}")
for seed, v in sorted(nearby.items()):
    print(f"  seed {seed}: E[nearby retained] at tau={TAU} is {v[TAU]:.0f}")
ref_seed = 23 if 23 in nearby else min(nearby)

print(f"\nTable 13 (seed {ref_seed}); T_key floor = two stages of Step 4 = 2^{T_STAGES:.2f}")
rows = threshold_table(LOG_N, LOG_P, BLOCK, K, TAUS, nearby[ref_seed], T_STAGES)
print_table(rows, "T_key")

succ = poisson_ge(2.0 ** (LOG_N + LOG_P), TAU)
lo = min(v[TAU] for v in nearby.values())
hi = max(v[TAU] for v in nearby.values())
fin = math.log2((1 + hi) * 2.0 ** (128 - K))
print(f"\ntau = {TAU}: Pr[c_* >= tau] = {succ:.5f}, success probability "
      f"{succ - (1 + hi) * 2.0 ** (128 - K - BLOCK):.4f}")
print(f"  E[N_wrong] = {lo:.0f} .. {hi:.0f} (= 2^{math.log2(hi):.1f}), "
      f"final stage 2^{fin:.1f}, T_key = 2^{math.log2(2**T_STAGES + 2**fin):.2f}")
print(f"  even 2^20 retained classes would give "
      f"T_key = 2^{math.log2(2**T_STAGES + 2.0**(20 + 128 - K)):.2f}")
print(f"  lower bound from the A_j/B_j construction: Pr[c_* >= {2*TAU-1}] = "
      f"{poisson_ge(32.0, 2*TAU-1):.4f}, E[N_wrong] >= {6*poisson_ge(32.0, 2*TAU-1):.3f}")
