#!/usr/bin/env python3
"""Derivation of the masks and numbers used in the key-recovery attack.

Self-contained (Python 3 standard library only).  Words are row-wise:
x = x{3}||x{2}||x{1}||x{0}, x{i} = bits 8i+7..8i, and column j of x is the
four-bit SSb input u_j = (x{3}_j, x{2}_j, x{1}_j, x{0}_j).

f(x) = P32(G0(x) ^ G1(x <<<_8 1)).  Column q of G0(x) ^ G1(x <<<_8 1) is
N0(u_q) ^ N1(u_{q-1}); P32 moves bit q of every word to position p8(q).
Since different columns use disjoint input bits, every statement below about
D_delta f is obtained by exhaustive enumeration over each 4-bit column, i.e.
it is exact (no sampling).

[1] SSb: 97 jointly realizable transitions vs. 157 admitted by the product of
    the marginal DDTs (Sect. V-A).
[2] Filters (Sect. V-C, Eqs. (5), (6)): bits of dx1 and dx0 that are fixed for
    every pair reaching (dx3, dy3) = (00900000, 008C0000).
[3] Key masks (Appendix): M^(1) for testing dx2 = 008C0000, M^(2) for the
    coordinates of x2 needed in the test of dx3, M_K0 = M^(1)|M^(2), and
    M_K1 for sk1.
[4] Numbers of Sect. V-C: lambda_right, surviving pairs, wrong classes
    reaching tau, success probability, time and memory complexities, and
    the comparison with the 15-round rectangle attack.
"""
import itertools
import math
import random

P8 = [7, 4, 1, 6, 3, 0, 5, 2]
P8INV = [P8.index(j) for j in range(8)]
M32 = 0xFFFFFFFF


def N0(u):                       # G0 on one column, u = (u3,u2,u1,u0)
    x3, x2, x1, x0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y0 = (x3 & x2) ^ x0
    y3 = (y0 & x1) ^ x3
    return (y3 << 3) | (x2 << 2) | (x1 << 1) | y0


def N1(u):                       # G1 on one column
    x3, x2, x1, x0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y2 = (x3 & x1) ^ x2
    y1 = (y2 & x0) ^ x1
    return (x3 << 3) | (y2 << 2) | (y1 << 1) | x0


def col(x, j):                   # u_j of a row-wise word
    return sum(((x >> (8 * i + j)) & 1) << i for i in range(4))


def from_cols(cols):             # inverse of col: cols[q] -> word (no P32)
    w = 0
    for q, u in enumerate(cols):
        for i in range(4):
            w |= ((u >> i) & 1) << (8 * i + q)
    return w


def p32(w):
    out = 0
    for i in range(4):
        for j in range(8):
            out |= ((w >> (8 * i + j)) & 1) << (8 * i + P8[j])
    return out


def f(x):                        # reference implementation of f
    cols = [col(x, j) for j in range(8)]
    z = [N0(cols[q]) ^ N1(cols[(q - 1) % 8]) for q in range(8)]
    return p32(from_cols(z))


def f_ref(x):                    # f from the word-level specification
    def word(v, i): return (v >> (8 * i)) & 0xFF
    def pack(a, b, c, d): return (a << 24) | (b << 16) | (c << 8) | d
    def G0(v):
        x3, x2, x1, x0 = (word(v, i) for i in (3, 2, 1, 0))
        y0 = (x3 & x2) ^ x0
        return pack((y0 & x1) ^ x3, x2, x1, y0)
    def G1(v):
        x3, x2, x1, x0 = (word(v, i) for i in (3, 2, 1, 0))
        y2 = (x3 & x1) ^ x2
        return pack(x3, y2, (y2 & x0) ^ x1, x0)
    rot = pack(*(((word(x, i) << 1) | (word(x, i) >> 7)) & 0xFF for i in (3, 2, 1, 0)))
    return p32(G0(x) ^ G1(rot))


def hexw(w):
    return f"{w:08X}"


random.seed(1)
assert all(f(x) == f_ref(x) for x in (random.getrandbits(32) for _ in range(20000)))
print("[0] column model of f agrees with the word-level specification (20000 inputs)")

# --------------------------------------------------------------------- [1]
joint = {(a, N0(u) ^ N0(u ^ a), N1(u) ^ N1(u ^ a)) for a in range(16) for u in range(16)}
m0 = {(a, N0(u) ^ N0(u ^ a)) for a in range(16) for u in range(16)}
m1 = {(a, N1(u) ^ N1(u ^ a)) for a in range(16) for u in range(16)}
marg = {(a, b, c) for (a, b) in m0 for (a2, c) in m1 if a == a2}
print(f"[1] SSb transitions: joint {len(joint)}, product of marginals {len(marg)}")


# --------------------------------------------------------------------- [2]
def deriv_fixed(delta):
    """Exact (mask, value) of the bits of D_delta f(x) constant over all x."""
    d = [col(delta, j) for j in range(8)]
    A0 = [{N0(u) ^ N0(u ^ d[j]) for u in range(16)} for j in range(8)]
    A1 = [{N1(u) ^ N1(u ^ d[j]) for u in range(16)} for j in range(8)]
    mask = val = 0
    for q in range(8):
        for i in range(4):
            b0 = {(v >> i) & 1 for v in A0[q]}
            b1 = {(v >> i) & 1 for v in A1[(q - 1) % 8]}
            if len(b0) == 1 and len(b1) == 1:
                pos = 8 * i + P8[q]
                mask |= 1 << pos
                val |= (b0.pop() ^ b1.pop()) << pos
    return mask, val


def deriv_values(delta):
    """Exact set of values of D_delta f(x) over all x (active columns only)."""
    d = [col(delta, j) for j in range(8)]
    need = [j for j in range(8) if d[j]]          # only active columns matter
    out = set()
    for us in itertools.product(range(16), repeat=len(need)):
        u = dict(zip(need, us))
        z = [0] * 8
        for j in need:
            z[j] ^= N0(u[j]) ^ N0(u[j] ^ d[j])            # G0 part of column j
            z[(j + 1) % 8] ^= N1(u[j]) ^ N1(u[j] ^ d[j])  # G1 part moves to j+1
        out.add(p32(from_cols(z)))
    return out


DX3, DY3 = 0x00900000, 0x008C0000
DX2 = DY3
S1 = {DX3 ^ v for v in deriv_values(DX2)}          # all reachable dx1
m1_, v1_ = deriv_fixed(DX2)
m1_mask, m1_val = m1_, (v1_ ^ DX3) & m1_
print(f"[2] dx1 = 00900000 ^ D_(008C0000) f(x2): {len(S1)} values,"
      f" fixed bits mask {hexw(m1_mask)} value {hexw(m1_val)}  (Eq. (6): B9FF37B9 / 001E0000)")
assert all(d & m1_mask == m1_val for d in S1)
# dx0 = dx2 ^ D_(dx1) f(x1) over all reachable dx1 and all x1
mask0, val0 = M32, None
for d1 in S1:
    mk, vl = deriv_fixed(d1)
    vl ^= DX2 & mk
    if val0 is None:
        mask0, val0 = mk, vl
    else:
        mask0 &= mk & ~(vl ^ val0)
        val0 &= mask0
val0 &= mask0
print(f"    dx0 = 008C0000 ^ D_(dx1) f(x1): fixed bits mask {hexw(mask0)} value {hexw(val0)}"
      f"  (Eq. (5): 81311081 / 00110000)")
print(f"    Hamming weights {bin(mask0).count('1')} and {bin(m1_mask).count('1')}")


# --------------------------------------------------------------------- [3]
def deriv_support(delta):
    """Input bits of x on which D_delta f(x) actually depends (exact)."""
    d = [col(delta, j) for j in range(8)]
    sup = 0
    for j in range(8):
        if not d[j]:
            continue
        for bit in range(4):
            dep = any((N0(u) ^ N0(u ^ d[j])) != (N0(u ^ (1 << bit)) ^ N0(u ^ (1 << bit) ^ d[j])) or
                      (N1(u) ^ N1(u ^ d[j])) != (N1(u ^ (1 << bit)) ^ N1(u ^ (1 << bit) ^ d[j]))
                      for u in range(16))
            if dep:
                sup |= 1 << (8 * bit + j)
    return sup


V1 = [m1_val | sum(((s >> k) & 1) << p for k, p in enumerate(
      [p for p in range(32) if not (m1_mask >> p) & 1])) for s in range(1 << 9)]
M1 = 0
for d1 in V1:                                    # every dx1 allowed by Eq. (6)
    M1 |= deriv_support(d1)
print(f"[3] M^(1) (sk0 bits for testing dx2 = 008C0000): {hexw(M1)}  (paper: DECE5EDE)")

# bits of x2 needed for the test of dx3 = dx1 ^ D_(008C0000) f(x2)
need_x2 = deriv_support(DX2)
print(f"    x2 bits needed for the test of dx3: {hexw(need_x2)}"
      f" -> M_K1 = {hexw(need_x2)}  (paper: 8C008C8C), weight {bin(need_x2).count('1')}")


def f_support_of_output_bits(outmask):
    """Input bits of x on which the output bits outmask of f(x) depend."""
    sup = 0
    for i in range(4):
        for k in range(8):
            if not (outmask >> (8 * i + k)) & 1:
                continue
            q = P8INV[k]                          # column of z feeding bit k
            for (c, fn) in ((q, N0), ((q - 1) % 8, N1)):
                for bit in range(4):
                    if any(((fn(u) >> i) & 1) != ((fn(u ^ (1 << bit)) >> i) & 1) for u in range(16)):
                        sup |= 1 << (8 * bit + c)
    return sup


M2 = f_support_of_output_bits(need_x2)          # x2 = x0 ^ f(x1) ^ sk1
print(f"    M^(2) (sk0 bits for the coordinates of x2): {hexw(M2)}  (paper: D9D9D9D9)")
MK0 = M1 | M2
print(f"    M_K0 = M^(1) | M^(2) = {hexw(MK0)}  (paper: DFDFDFDF), weight {bin(MK0).count('1')};"
      f" total {bin(MK0).count('1') + bin(need_x2).count('1')} guessed bits")

# --------------------------------------------------------------------- [4]
def pois_sf(k, lam):              # Pr[Pois(lam) > k], summed over the tail
    return sum(math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1)) for i in range(k + 1, k + 400))


p_diff = 2.0**-57
lam_r = 2**62 * p_diff
surv = 2**62 * 2.0**-8 * 2.0**-23
lam_w = surv * 2.0**-33
wrong = (2**37 - 1) * pois_sf(15, lam_w)
print(f"[4] lambda_right = 2^62 * 2^-57 = {lam_r:g};  surviving pairs 2^{math.log2(surv):g};"
      f"  lambda_wrong = 2^{math.log2(lam_w):g}")
print(f"    Pr[Pois(32) >= 16] = {pois_sf(15, 32.0):.5f};"
      f"  (2^37-1) Pr[Pois(2^-2) >= 16] = 2^{math.log2(wrong):.1f}")
T = 2 * 2 / 16 * 2**62 + 2 / 16 * 2**54 + 2 * 2 / 16 * 2**31 * 2**37 + 2**91
print(f"    T = 2^{math.log2(T):.4f};  worst case 27 * 2^91 = 2^{math.log2(27 * 2**91):.2f}")
M = 2**37 + 2 * 2**31 * 64 / 8
print(f"    M = 2^37 + 2 * 2^31 * 64/8 = 2^{math.log2(M):.2f} bytes")
s = 2
Pmin = s * 2.0**-61
print(f"    rectangle: s = N^2 2^-65 P with N = 2^55, P = 2^-43.23 -> s = {2**110 * 2**-65 * 2**-43.23:.2f};"
      f"  N <= 2^63 -> P >= s 2^-61 = 2^{math.log2(Pmin):g};"
      f"  margin 2^{math.log2(Pmin) + 43.23:.2f}")

ok = (len(joint), len(marg)) == (97, 157) and (m1_mask, m1_val) == (0xB9FF37B9, 0x001E0000) \
    and (mask0, val0) == (0x81311081, 0x00110000) and M1 == 0xDECE5EDE and M2 == 0xD9D9D9D9 \
    and MK0 == 0xDFDFDFDF and need_x2 == 0x8C008C8C
print("\nALL VALUES AGREE WITH THE PAPER" if ok else "\nMISMATCH WITH THE PAPER")
