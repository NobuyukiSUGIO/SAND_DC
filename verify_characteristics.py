#!/usr/bin/env python3
"""Independent verification of the published SAND-64 differential characteristics.

This script does NOT use the SAT model.  For every round transition

    (dX^r, dY^r) -> (dX^{r+1}, dY^{r+1})

it counts the exact number of 32-bit round inputs that realize the transition,
using a transfer-matrix product over the eight SSb positions (see
verify_sand64_key_recovery.round_transition_count).  The weight of the
transition is 32 - log2(count); a count of zero proves the transition to be
impossible.  Because the transfer matrix walks the eight synthetic S-boxes
SSb = (G0, G1) jointly, it enforces the constraint that G0 and G1 consume the
same nibble -- exactly the constraint that a model built from two independent
DDTs would miss.

Usage
    python3 verify_characteristics.py                       # verify everything shipped in results/
    python3 verify_characteristics.py FILE EXPECTED_WEIGHT  # verify one file
"""
import importlib.util
import math
import multiprocessing as mp
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import trail_io

_spec = importlib.util.spec_from_file_location(
    "sand64_verify", os.path.join(HERE, "verify_sand64_key_recovery.py"))
V = importlib.util.module_from_spec(_spec)
sys.modules["sand64_verify"] = V
_spec.loader.exec_module(V)


def exact_round_weights(trail):
    """Exact weight of every round transition, or None if some transition is impossible."""
    weights = []
    for r in range(len(trail) - 1):
        count = V.round_transition_count(*trail[r], *trail[r + 1])
        if count == 0:
            return None
        weights.append(32.0 - math.log2(count))
    return weights


def verify_file(path, expected_weight, expected_count=None, workers=None):
    name = os.path.basename(path)
    trails = trail_io.read_trails(path)
    if not trails:
        print(f"[{name}] no characteristic found -- FAIL")
        return False

    rounds = {len(t) - 1 for t in trails}
    print(f"[{name}] {len(trails)} characteristic(s), {sorted(rounds)} round(s)")

    ok = True
    if len(rounds) != 1:
        print(f"   inconsistent round numbers {sorted(rounds)} -- FAIL")
        ok = False

    distinct = len({tuple(t) for t in trails})
    if distinct != len(trails):
        print(f"   {len(trails) - distinct} duplicate record(s) -- FAIL")
        ok = False

    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    if len(trails) > 1:
        with mp.Pool(processes=workers) as pool:
            res = pool.map(exact_round_weights, trails, chunksize=8)
    else:
        res = [exact_round_weights(trails[0])]

    impossible = [i for i, r in enumerate(res) if r is None]
    if impossible:
        print(f"   {len(impossible)} characteristic(s) contain an impossible round "
              f"transition (first: #{impossible[0] + 1}) -- FAIL")
        ok = False

    sums = Counter(round(sum(r), 6) for r in res if r is not None)
    for w in sorted(sums):
        mark = "OK" if w == expected_weight else "FAIL (expected 2^-%g)" % expected_weight
        print(f"   exact probability 2^-{w:g} : {sums[w]} characteristic(s) -- {mark}")
        if w != expected_weight:
            ok = False

    if expected_count is not None and len(trails) != expected_count:
        print(f"   expected {expected_count} characteristics, found {len(trails)} -- FAIL")
        ok = False

    if len(trails) > 1:
        ins = {(t[0][0], t[0][1]) for t in trails}
        outs = {(t[-1][0], t[-1][1]) for t in trails}
        pairs = Counter(((t[0][0], t[0][1]), (t[-1][0], t[-1][1])) for t in trails)
        print(f"   distinct input differences  : {len(ins)}")
        print(f"   distinct output differences : {len(outs)}")
        print(f"   distinct (dP, dC) pairs     : {len(pairs)}")
        mult = Counter(pairs.values())
        for k in sorted(mult):
            print(f"     {mult[k]} pair(s) realized by {k} characteristic(s)"
                  f"  => differential probability >= 2^-{expected_weight - math.log2(k):g}")
        vecs = Counter(tuple(int(round(x)) for x in r) for r in res if r is not None)
        print("   round-weight vectors:")
        for v, c in sorted(vecs.items(), key=lambda kv: -kv[1]):
            print(f"     {list(v)}  sum={sum(v)} : {c}")

    return ok


def main():
    if len(sys.argv) >= 3:
        ok = verify_file(sys.argv[1], float(sys.argv[2]))
    else:
        res_dir = os.path.join(HERE, "results")
        jobs = [(os.path.join(res_dir, "trails_w58_all.txt"), 58.0, 12288),
                (os.path.join(res_dir, "characteristic_13r_w64.txt"), 64.0, 1)]
        ok = True
        for path, w, count in jobs:
            if not os.path.exists(path):
                print(f"[{os.path.basename(path)}] missing -- skipped")
                continue
            ok &= verify_file(path, w, count)
            print()
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
