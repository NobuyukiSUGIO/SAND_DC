#!/usr/bin/env python3
"""Exact transfer-matrix computation of SAND-64 round-transition counts.

round_transition_count(dx, dy, dx_next, dy_next) returns the exact number of
32-bit round inputs x realizing the round transition
(dx, dy) -> (dx_next, dy_next); the weight is 32 - log2(count), and a count
of zero proves the transition to be impossible.  The computation multiplies
eight 16x16 transfer matrices induced by the joint DDT of the synthetic S-box
SSb = (G0, G1) along the cyclic coupling of the eight SSb positions, so the
constraint that G0 and G1 consume the same nibble is enforced exactly.  It
does not use the SAT model.

The state words are row-wise (x{3}||x{2}||x{1}||x{0}); they are transposed to
the nibble-wise SSb representation used by SAND.
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
