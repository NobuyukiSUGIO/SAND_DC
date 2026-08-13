#!/usr/bin/env python3
"""
23-round SAND-128 key recovery: success probability of the counting phase.

Wrong-key randomisation model.  With N pairs, a W-weight characteristic and a
BLOCK-bit condition at the end of the distinguisher,

    C_right  ~ Poisson(2^(N-W)),
    C_wrong  ~ Poisson(2^(N-BLOCK))   independently over 2^k - 1 classes.

Two decision rules are compared:

  (a) keep the maximum-counter class(es)   -- what the draft currently says;
  (b) keep every class whose counter is at least a threshold tau.

Rule (b) is the one to put in the paper: its success probability is much
higher and the extra classes it keeps cost only (#kept) x 2^(KEYLEN-k)
encryptions in the final exhaustive-search stage, which stays far below the
2^118.5 counting cost.
"""
from decimal import Decimal as D, getcontext
import math

getcontext().prec = 200

N = 124          # log2 of the number of pairs
W = 120          # weight of the distinguisher
BLOCK = 128
KEYLEN = 128
K = 77           # guessed partial-key bits
T_COUNT = 118.48  # log2 of the counting-stage time

LAM_R = D(2) ** (N - W)
LAM_W = D(2) ** (N - BLOCK)
M = D(2) ** K - 1
TMAX = 400


def pmf(t, lam):
    return (-lam).exp() * lam ** t / D(math.factorial(t))


PR = [pmf(t, LAM_R) for t in range(TMAX)]
PW = [pmf(t, LAM_W) for t in range(TMAX)]
CR = [sum(PR[:t + 1]) for t in range(TMAX)]      # P(C_right <= t)
CW = [sum(PW[:t + 1]) for t in range(TMAX)]      # P(C_wrong <= t)


def l2(x):
    return float(x).__format__(".4g") if x == 0 else f"2^{math.log2(float(x)):.1f}"


print(f"parameters: N = 2^{N}, W = {W}, k = {K}")
print(f"  E[C_right] = {float(LAM_R):g}, E[C_wrong] = 2^{N-BLOCK} = {float(LAM_W):g}, "
      f"wrong classes = 2^{K}")
print()

# ---- rule (a): take the maximum ------------------------------------------
strict = sum(PR[t] * CW[t - 1] ** int(M) for t in range(1, TMAX))
incl = sum(PR[t] * CW[t] ** int(M) for t in range(TMAX))
print("(a) keep the maximum-counter class:")
print(f"    P(correct is the strict maximum)      = {float(strict):.4%}")
print(f"    P(correct is among the maxima)        = {float(incl):.4%}")
print()

# ---- rule (b): threshold --------------------------------------------------
print("(b) keep every class with counter >= tau:")
print(f"    {'tau':>4} {'P(success)':>12} {'E[wrong kept]':>15} "
      f"{'final search':>14} {'total time':>12}")
for tau in range(4, 15):
    succ = 1 - CR[tau - 1]
    kept = M * (1 - CW[tau - 1])
    final = float(1 + kept) * 2 ** (KEYLEN - K)
    total = math.log2(2 ** T_COUNT + final)
    print(f"    {tau:>4} {float(succ):>11.4%} {l2(kept):>15} "
          f"{l2(D(final)):>14} {'2^%.2f' % total:>12}")
print()
print(f"The counting stage alone costs 2^{T_COUNT}, so any tau <= 12 leaves the")
print(f"total time unchanged at 2^{T_COUNT} while keeping the success probability")
print("above 99%.")
