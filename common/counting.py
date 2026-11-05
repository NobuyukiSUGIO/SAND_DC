"""
Statistical model of the counting phase (Section 7.3.1 of the paper).

  c_*      ~ Poisson(N p)                 correct class
  c_kappa  ~ Poisson(N 2^-b)              background wrong class
  c_kappa  ~ Poisson(N p rho + N 2^-b)    nearby wrong class, measured by
                                          nearby/nearby_classes.py

Every retained class is searched exhaustively, so

  Pr[success]   = Pr[c_* >= tau] - Pr[false positive in the final check]
  final stage   = (1 + E[N_wrong]) * 2^(128 - k)
  E[N_wrong]    = E[background retained] + E[nearby retained].

Standard library only.
"""
import glob
import json
import math
import os

NEARBY_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "nearby",
                          "results")


def poisson_ge(lam, k):
    """Pr[Poisson(lam) >= k], accurate also in the far tail."""
    if lam < k:
        t = math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))
        s, i = 0.0, k
        while t > 1e-30 * max(s, 1e-300):
            s += t
            i += 1
            t *= lam / i
        return s
    s, t = 0.0, math.exp(-lam)
    for i in range(k):
        s += t
        t *= lam / (i + 1)
    return max(0.0, 1.0 - s)


def load_nearby(variant):
    """{seed: {tau: E[nearby retained]}} from nearby/results/*.json."""
    out = {}
    for path in sorted(glob.glob(os.path.join(NEARBY_DIR, f"nearby_sand{variant}_seed*.json"))):
        d = json.load(open(path))
        out[d["seed"]] = {int(t): v for t, v in d["nearby"].items()}
    return out


def threshold_table(log_n, log_p, block, k, taus, nearby, log_t_floor, keylen=128):
    """Rows (tau, Pr[c_* >= tau], E[bg], E[nearby], log2 final, log2 total)."""
    lam_r, lam_b = 2.0 ** (log_n + log_p), 2.0 ** (log_n - block)
    rows = []
    for tau in taus:
        succ = poisson_ge(lam_r, tau)
        bg = (2.0 ** k - 1) * poisson_ge(lam_b, tau)
        nb = nearby.get(tau, float("nan"))
        final = (1 + bg + nb) * 2.0 ** (keylen - k)
        total = math.log2(2.0 ** log_t_floor + final)
        rows.append((tau, succ, bg, nb, math.log2(final), total))
    return rows


def print_table(rows, floor_name):
    print(f"  {'tau':>4} {'Pr[c*>=tau]':>12} {'E[bg]':>11} {'E[nearby]':>10} "
          f"{'final':>8} {floor_name:>8}")
    for tau, succ, bg, nb, fin, tot in rows:
        bgs = f"2^{math.log2(bg):.1f}" if bg > 0 else "0"
        print(f"  {tau:>4} {succ:>12.5f} {bgs:>11} {nb:>10.1f} "
              f"{'2^%.1f' % fin:>8} {'2^%.1f' % tot:>8}")
