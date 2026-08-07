#!/usr/bin/env python3
"""Independent checks for the SAND-64 differential/key-recovery draft.

Checks performed:
1. Build the joint 4x8 synthetic S-box SSb=(G0,G1) from the SAND spec.
2. Recompute the exact probability of every round transition in Table III.
3. Evaluate hypothetical key-ranking success under a corrected valid
   characteristic, using the wrong-key randomization/Poisson approximation.

The state words in the draft are row-wise.  The exact differential calculation
transposes them to the nibble-wise SSb representation used by SAND.
"""
from __future__ import annotations

import math
from typing import List, Sequence, Tuple


def g0(x: int) -> int:
    x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
    y0 = (x3 & x2) ^ x0
    y1 = x1
    y2 = x2
    y3 = (y0 & x1) ^ x3
    return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0


def g1(x: int) -> int:
    x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
    y3 = x3
    y2 = (x3 & x1) ^ x2
    y1 = (y2 & x0) ^ x1
    y0 = x0
    return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0


def build_joint_ddt() -> List[List[int]]:
    ddt = [[0] * 256 for _ in range(16)]
    for din in range(16):
        for x in range(16):
            dg0 = g0(x) ^ g0(x ^ din)
            dg1 = g1(x) ^ g1(x ^ din)
            ddt[din][(dg0 << 4) | dg1] += 1
    return ddt


DDT = build_joint_ddt()
P8 = [7, 4, 1, 6, 3, 0, 5, 2]
PERM = [4 * P8[i // 4] + (i % 4) for i in range(32)]


def row_to_ssb(word: int) -> int:
    """Transpose x{3}||x{2}||x{1}||x{0} into eight contiguous nibbles."""
    rows = [word & 0xFF, (word >> 8) & 0xFF, (word >> 16) & 0xFF, (word >> 24) & 0xFF]
    out = 0
    for col in range(8):
        nib = (
            (((rows[3] >> col) & 1) << 3)
            | (((rows[2] >> col) & 1) << 2)
            | (((rows[1] >> col) & 1) << 1)
            | ((rows[0] >> col) & 1)
        )
        out |= nib << (4 * col)
    return out


def inverse_pn(diff: int) -> int:
    """Return the pre-permutation XOR-compression difference."""
    out = 0
    # The official model uses permG[PERM[i]] = xorG[i].
    for i in range(32):
        out |= ((diff >> PERM[i]) & 1) << i
    return out


def nibbles(word: int) -> List[int]:
    return [(word >> (4 * i)) & 0xF for i in range(8)]


def matmul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> List[List[int]]:
    c = [[0] * 16 for _ in range(16)]
    for i in range(16):
        for k in range(16):
            aik = a[i][k]
            if not aik:
                continue
            for j in range(16):
                bkj = b[k][j]
                if bkj:
                    c[i][j] += aik * bkj
    return c


def round_transition_count(
    dx: int, dy: int, dx_next: int, dy_next: int
) -> int:
    """Number of 32-bit x values realizing one SAND round transition.

    The denominator is 2^32.  A zero count proves that the transition is
    impossible.  The transfer-matrix trace enforces the joint G0/G1 output
    constraint for every shared 4-bit input nibble.
    """
    dx = row_to_ssb(dx)
    dy = row_to_ssb(dy)
    dx_next = row_to_ssb(dx_next)
    dy_next = row_to_ssb(dy_next)
    if dy_next != dx:
        return 0

    target_xor = inverse_pn(dy ^ dx_next)
    din = nibbles(dx)
    z = nibbles(target_xor)

    matrices: List[List[List[int]]] = []
    for i in range(8):
        # rotl(G1,4): z_i = dG0_i xor dG1_{i-1}.
        m = [[0] * 16 for _ in range(16)]
        for dout, count in enumerate(DDT[din[i]]):
            if not count:
                continue
            dg0, dg1 = (dout >> 4) & 0xF, dout & 0xF
            previous_dg1 = z[i] ^ dg0
            m[previous_dg1][dg1] += count
        matrices.append(m)

    product = [[int(i == j) for j in range(16)] for i in range(16)]
    for m in matrices:
        product = matmul(product, m)
    return sum(product[i][i] for i in range(16))


TRAIL: List[Tuple[int, int]] = [
    (0x00460000, 0x00F00000),
    (0x00840000, 0x00460000),
    (0x00800000, 0x00840000),
    (0x00000000, 0x00800000),
    (0x00800000, 0x00000000),
    (0x00840000, 0x00800000),
    (0x00460000, 0x00840000),
    (0x00F04000, 0x00460000),
    (0x00CA0400, 0x00F04000),
    (0x000A8000, 0x00CA0400),
    (0x00100000, 0x000A8000),
    (0x00038000, 0x00100000),
    (0x04929600, 0x00038000),
]
CLAIMED_WEIGHTS = [6, 4, 2, 0, 2, 4, 6, 9, 10, 7, 2, 6]


def verify_trail() -> None:
    print("Exact joint-SSb verification of Table III")
    print("round  claimed  recomputed")
    print("-----  -------  ----------")
    for r in range(12):
        count = round_transition_count(*TRAIL[r], *TRAIL[r + 1])
        if count == 0:
            result = "IMPOSSIBLE"
        else:
            weight = 32.0 - math.log2(count)
            result = f"{weight:g}"
        print(f"{r:>2}->{r+1:<2}  {CLAIMED_WEIGHTS[r]:>7}  {result:>10}")

    # A compact certificate for why the separate-marginal-DDT model can accept
    # the claimed weight-9 transition although the exact joint model rejects it.
    print("\nRound 7->8 marginal-model conflict certificate:")
    print("  one minimum-weight assignment admitted by separate DDTs uses")
    print("  input nibble difference = 0x6 and (dG0,dG1) = (0x6,0x2)")
    print(f"  but the joint entry DDT_SSb[0x6][0x62] = {DDT[0x6][0x62]}")
    print("  exact 32-bit round-transition count = 0")


def poisson_rank_stats(lambda_wrong: float, lambda_right: float, candidates: int) -> Tuple[float, float, float]:
    """Strict-win, max-list inclusion, and random-tie top-1 success.

    Requires mpmath for adequate precision with 2^52 candidates.
    """
    try:
        import mpmath as mp
    except ImportError as exc:
        raise SystemExit("Install mpmath to run the statistical section: pip install mpmath") from exc

    mp.mp.dps = 100
    m = mp.mpf(candidates - 1)
    lw = mp.mpf(str(lambda_wrong))
    lr = mp.mpf(str(lambda_right))

    def pmf(k: int, lam):
        return mp.e ** (-lam) * lam**k / mp.factorial(k)

    def cdf(k: int, lam):
        if k < 0:
            return mp.mpf("0")
        return mp.fsum(pmf(i, lam) for i in range(k + 1))

    strict = mp.mpf("0")
    inclusive = mp.mpf("0")
    random_tie = mp.mpf("0")
    for t in range(100):
        pr = pmf(t, lr)
        below = cdf(t - 1, lw)
        at_most = cdf(t, lw)
        equal = pmf(t, lw)
        strict += pr * below**m
        inclusive += pr * at_most**m
        random_tie += pr * (
            (at_most ** (m + 1) - below ** (m + 1)) / ((m + 1) * equal)
        )
    return float(strict), float(inclusive), float(random_tie)


def conditional_statistics() -> None:
    n_pairs = 2**62
    key_candidates = 2**52
    p_right = 2.0**-58
    lambda_right = n_pairs * p_right

    print("\nHypothetical statistics after replacing Table III by a VALID p=2^-58 trail")
    for tested_bits in (64, 68):
        p_wrong = 2.0 ** (-tested_bits)
        lambda_wrong = n_pairs * p_wrong
        strict, inclusive, tie = poisson_rank_stats(lambda_wrong, lambda_right, key_candidates)
        print(f"\nWrong-key model: {tested_bits} independent tested bits")
        print(f"  E[C_right] = {lambda_right:g}")
        print(f"  E[C_wrong] = {lambda_wrong:g}")
        print(f"  probability-ratio S/N = 2^{tested_bits-58} = {p_right/p_wrong:g}")
        print(f"  P(correct strictly highest) = {strict:.9%}")
        print(f"  P(correct in maximum-counter list) = {inclusive:.9%}")
        print(f"  P(top-1 with random tie-breaking) = {tie:.9%}")


if __name__ == "__main__":
    verify_trail()
    conditional_statistics()
