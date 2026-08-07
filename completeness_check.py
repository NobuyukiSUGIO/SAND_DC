#!/usr/bin/env python3
"""Proof that the published list of weight-58 characteristics is complete.

The 12-round SAT model is built with the cardinality constraint "total weight
is exactly 58" and NO partitioning constraint of any kind.  Every characteristic
of results/trails_w58_all.txt is then added as a blocking clause over the state
differences (dY^0, dX^0, dX^1, ..., dX^12), which together determine a
characteristic uniquely (dY^r = dX^{r-1} holds in the forward model).

If the resulting formula is UNSATISFIABLE, no 12-round characteristic of weight
58 exists outside the published list.  This argument is independent of how the
enumeration itself was organised.

Usage
    python3 completeness_check.py [TRAIL_FILE]
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import trail_io
from enumerate_w58 import key_vars, run_kissat, extract_trail, round_weights
import enum_by_wvec as W

NUM_ROUNDS = 12
DEFAULT_TRAILS = os.path.join(HERE, "results", "trails_w58_all.txt")
TMP = os.environ.get("SAND_TMPDIR", "/dev/shm")


def blocking_clause(trail, kv):
    """Forbid exactly this characteristic (negation of its state-difference assignment)."""
    bits = [(trail[0][1] >> (31 - i)) & 1 for i in range(32)]          # dY^0
    for r in range(NUM_ROUNDS + 1):
        bits += [(trail[r][0] >> (31 - i)) & 1 for i in range(32)]     # dX^0 .. dX^12
    return [(-v if b else v) for v, b in zip(kv, bits)]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_TRAILS
    trails = trail_io.read_trails(path, NUM_ROUNDS)
    print(f"loaded {len(trails)} characteristics from {os.path.basename(path)}", flush=True)
    if len({tuple(t) for t in trails}) != len(trails):
        sys.exit("the trail file contains duplicates")

    model, variables = W.build(with_counters=False)   # weight == 58, no partitioning
    kv = key_vars(variables, NUM_ROUNDS)
    for t in trails:
        model.clauses.append(blocking_clause(t, kv))

    num_vars = max(max(abs(l) for l in c) for c in model.clauses)
    print(f"CNF: {num_vars} variables, {len(model.clauses)} clauses", flush=True)

    cnf = os.path.join(TMP, "sand64_complete.cnf")
    with open(cnf, "w") as f:
        f.write(f"p cnf {num_vars} {len(model.clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in model.clauses))

    t0 = time.time()
    sat, pos = run_kissat(model.kissat_path, cnf)
    dt = time.time() - t0
    os.remove(cnf)

    if sat:
        tr = extract_trail(pos, variables, NUM_ROUNDS)
        print(f"*** SATISFIABLE ({dt:.0f}s): the list is INCOMPLETE. Missing characteristic:")
        for r, (dx, dy) in enumerate(tr):
            print(f"  {r:>3} {trail_io.fmt(dx)} {trail_io.fmt(dy)}")
        print("  round weights:", round_weights(pos, variables, NUM_ROUNDS))
        sys.exit(1)

    print(f"UNSATISFIABLE ({dt:.0f}s)")
    print(f"=> completeness proved: exactly {len(trails)} 12-round differential "
          f"characteristics of SAND-64 have weight 58 (p = 2^-58); no other exists.")


if __name__ == "__main__":
    main()
