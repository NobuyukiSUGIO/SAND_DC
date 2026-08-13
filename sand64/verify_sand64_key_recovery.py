#!/usr/bin/env python3
"""
Independent verification of the 17-round key-recovery attack on SAND-64
described in paper/SAND_DC_IEICE/main.tex.

SAND-64 is re-implemented from the specification (four 8-bit rows, p_8, the
row-wise rotation x <<<_{n/4} 1); nothing here is imported from the author's
search code.  The companion script for the other variant is
verify_sand128_key_recovery.py.

Usage:  python3 verify_sand64_key_recovery.py
"""
import collections
import itertools
import json
import math
import os
import random
import sys

random.seed(20260811)

P8 = [7, 4, 1, 6, 3, 0, 5, 2]
M32 = (1 << 32) - 1
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

ok_all = True


def check(name, cond, extra=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK ' if cond else 'FAIL'}] {name}{(' : ' + extra) if extra else ''}")
    return cond


# ---------------------------------------------------------------- primitives
def G0(u):
    u3, u2, u1, u0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y0 = (u3 & u2) ^ u0
    y3 = (y0 & u1) ^ u3
    return (y3 << 3) | (u2 << 2) | (u1 << 1) | y0


def G1(u):
    u3, u2, u1, u0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y2 = (u3 & u1) ^ u2
    y1 = (y2 & u0) ^ u1
    return (u3 << 3) | (y2 << 2) | (y1 << 1) | u0


DDT = [[0] * 256 for _ in range(16)]
for _d in range(16):
    for _u in range(16):
        DDT[_d][((G0(_u) ^ G0(_u ^ _d)) << 4) | (G1(_u) ^ G1(_u ^ _d))] += 1

PTAB = [0] * 256
for _v in range(256):
    _o = 0
    for _j in range(8):                       # input column j -> column 7-p8[7-j]
        if (_v >> (7 - _j)) & 1:
            _o |= 1 << P8[7 - _j]
    PTAB[_v] = _o


def F(x):
    """f(x) = P_n(G0(x) XOR G1(x <<<_{n/4} 1)), row-parallel over the four rows."""
    r3, r2, r1, r0 = (x >> 24) & 0xFF, (x >> 16) & 0xFF, (x >> 8) & 0xFF, x & 0xFF
    a0 = (r3 & r2) ^ r0
    a3 = (a0 & r1) ^ r3
    s3 = ((r3 << 1) | (r3 >> 7)) & 0xFF       # each row rotated left by one
    s2 = ((r2 << 1) | (r2 >> 7)) & 0xFF
    s1 = ((r1 << 1) | (r1 >> 7)) & 0xFF
    s0 = ((r0 << 1) | (r0 >> 7)) & 0xFF
    b2 = (s3 & s1) ^ s2
    b1 = (b2 & s0) ^ s1
    return (PTAB[a3 ^ s3] << 24) | (PTAB[r2 ^ b2] << 16) | \
           (PTAB[r1 ^ b1] << 8) | PTAB[a0 ^ s0]


def enc(x, y, sk):                            # ordinary round, with the swap
    return (y ^ F(x) ^ sk) & M32, x


def dec(xn, yn, sk):
    x = yn
    return x & M32, (xn ^ F(x) ^ sk) & M32


def enc_last(x, y, sk):                       # last round, swap omitted
    return x, (y ^ F(x) ^ sk) & M32


def dec_last(x, y, sk):
    return x, (y ^ F(x) ^ sk) & M32


def cols(w):
    return [((w >> (31 - j)) & 1) << 3 | ((w >> (31 - (j + 8))) & 1) << 2 |
            ((w >> (31 - (j + 16))) & 1) << 1 | ((w >> (31 - (j + 24))) & 1)
            for j in range(8)]


def uncols(c):
    w = 0
    for j in range(8):
        for k, sh in enumerate((3, 2, 1, 0)):
            if (c[j] >> sh) & 1:
                w |= 1 << (31 - (j + 8 * k))
    return w


def p_inv(w):
    c = cols(w)
    return uncols([c[7 - P8[7 - j]] for j in range(8)])


def exact_count(dx_in, dy_in, dx_out):
    """number of 32-bit round inputs realizing the transition, by a transfer
    matrix over the eight SSb positions (joint N0/N1 constraint enforced)"""
    T = cols(p_inv(dx_out ^ dy_in))
    din = cols(dx_in)
    total = 0
    for d1_first in range(16):
        row = DDT[din[0]]
        cur = {}
        for d0 in range(16):
            c = row[(d0 << 4) | d1_first]
            if c:
                cur[d0] = cur.get(d0, 0) + c
        for m in range(1, 8):
            nxt = {}
            row = DDT[din[m]]
            for d0_prev, ways in cur.items():
                need = T[m - 1] ^ d0_prev
                for d0 in range(16):
                    c = row[(d0 << 4) | need]
                    if c:
                        nxt[d0] = nxt.get(d0, 0) + ways * c
            cur = nxt
            if not cur:
                break
        for d0_last, ways in cur.items():
            if T[7] ^ d0_last == d1_first:
                total += ways
    return total


rnd = lambda: random.getrandbits(32)

# ------------------------------------------------------- paper's constants
TRAIL = [(0x00900000, 0x008C0000), (0x00010000, 0x00900000), (0x00000000, 0x00010000),
         (0x00010000, 0x00000000), (0x00900000, 0x00010000), (0x008C0100, 0x00900000),
         (0x000E1000, 0x008C0100), (0x00950000, 0x000E1000), (0x00510000, 0x00950000),
         (0x00280000, 0x00510000), (0x00380000, 0x00280000), (0x00480000, 0x00380000),
         (0x00540000, 0x00480000)]
W = [4, 2, 0, 2, 4, 9, 9, 8, 6, 4, 6, 4]          # per-round weights, total 58
A, B = 0x00540000, 0x00480000                     # (dX^16, dY^16), ciphertext side
DX4, DY4 = 0x00900000, 0x008C0000                 # (dX^4 , dY^4 ), plaintext side
FM, FV = 0x81311081, 0x00110000                   # key-free filter on dX^1, 8 bits
S1M, S1V = 0xB9FF37B9, 0x001E0000                 # Stage 1 test on dX^2, 23 bits
MK0, MK1, MK2 = 0xFFFFFFFF, 0xDFDFDFDF, 0x8C008C8C
MK1a, MK1b = 0xDECE5EDE, 0xD9D9D9D9

print("=" * 78)
print("0.  self-consistency of the re-implementation")
print("=" * 78)
check("P_8 is a bit permutation",
      sorted(p_inv(1 << i) for i in range(32)) == sorted(1 << i for i in range(32)))
check("ordinary round is invertible",
      all(dec(*enc(x := rnd(), y := rnd(), k := rnd()), k) == (x, y) for _ in range(300)))
check("last round (swap omitted) is invertible",
      all(dec_last(*enc_last(x := rnd(), y := rnd(), k := rnd()), k) == (x, y)
          for _ in range(300)))
joint = sum(1 for a in range(16) for o in range(256) if DDT[a][o])
m0 = [[0] * 16 for _ in range(16)]
m1 = [[0] * 16 for _ in range(16)]
for a in range(16):
    for u in range(16):
        m0[a][G0(u) ^ G0(u ^ a)] += 1
        m1[a][G1(u) ^ G1(u ^ a)] += 1
marg = sum(1 for a in range(16) for b in range(16) for c in range(16) if m0[a][b] and m1[a][c])
check("joint SSb-DDT has 97 triples, product of marginals 157",
      joint == 97 and marg == 157, f"{joint} / {marg}")
check("DDT_SSb[0x6][0x62] = 0 (paper's example)", DDT[0x6][0x62] == 0)

print()
print("=" * 78)
print("1.  the 12-round characteristic (Table 4)")
print("=" * 78)
print(f"{'r':>3} {'dX^r':>10} {'dY^r':>10} {'paper':>6} {'exact':>6} {'count':>12}")
tot, link, wok = 0, True, True
for r in range(1, 13):
    dxp, dyp = TRAIL[r - 1]
    dxc, dyc = TRAIL[r]
    link &= (dyc == dxp)
    cnt = exact_count(dxp, dyp, dxc)
    w = 32 - cnt.bit_length() + 1 if cnt else None
    wok &= (w == W[r - 1]) and cnt > 0 and (cnt & (cnt - 1)) == 0
    tot += w if w else 0
    print(f"{r:>3} {dxc:08X} {dyc:08X} {W[r-1]:>6} {str(w):>6} {cnt:>12}")
check("Feistel linkage dY^r = dX^(r-1) holds for every round", link)
check("no impossible transition and every count is a power of two", wok)
check("exact total weight is 58, i.e. p = 2^-58", tot == 58, str(tot))

print("\n  Monte-Carlo check of each round transition (2^18 random inputs):")
N = 1 << 18
mc = True
for r in range(1, 13):
    dxp, dyp = TRAIL[r - 1]
    tgt = TRAIL[r][0] ^ dyp
    hits = 0
    for _ in range(N):
        x = rnd()
        if F(x) ^ F(x ^ dxp) == tgt:
            hits += 1
    est = math.log2(hits / N) if hits else None
    good = est is not None and abs(est + W[r - 1]) < 0.5
    mc &= good
    print(f"    round {r:>2}: claimed 2^-{W[r-1]:<2} measured "
          f"{('2^%.2f' % est) if est is not None else 'no hit':>9}")
check("all measured single-round DPs agree with the claimed weights", mc)

print("\n  composition over several rounds, with random round keys:")
for lo, hi in ((1, 5), (0, 5)):
    w = sum(W[lo:hi])
    n = 1 << (w + 11)
    hits = 0
    for _ in range(n):
        x, y = rnd(), rnd()
        xs, ys = x ^ TRAIL[lo][0], y ^ TRAIL[lo][1]
        for t in range(lo, hi):
            k = rnd()
            x, y = enc(x, y, k)
            xs, ys = enc(xs, ys, k)
            if (x ^ xs, y ^ ys) != TRAIL[t + 1]:
                break
        else:
            hits += 1
    est = math.log2(hits / n) if hits else None
    print(f"    rounds {lo}..{hi} (weight {w}): measured "
          f"{('2^%.2f' % est) if est is not None else 'none'}   expected 2^-{w}")
    check(f"segment {lo}..{hi} follows the whole path with the predicted probability",
          est is not None and abs(est + w) < 0.5)

print()
print("=" * 78)
print("2.  the enumeration of the optimal characteristics (Tables 5, 6)")
print("=" * 78)
path = os.path.join(RESULTS, "trails_w58_all.txt")
if not os.path.exists(path):
    print(f"  [skipped] {path} not found")
else:
    trails, cur = [], None
    for line in open(path):
        t = line.strip()
        if t.startswith("#"):
            if cur:
                trails.append(cur)
            cur = [] if t.startswith("# characteristic") else None
        elif t and cur is not None:
            p = t.split()
            cur.append((int(p[1] + p[2], 16), int(p[3] + p[4], 16)))
    if cur:
        trails.append(cur)
    trails = [t for t in trails if len(t) == 13]
    check("12288 characteristics, pairwise distinct",
          len(trails) == 12288 and len(set(map(tuple, trails))) == 12288, str(len(trails)))
    ins = {(t[0][0], t[0][1]) for t in trails}
    outs = {(t[12][0], t[12][1]) for t in trails}
    pairs = collections.Counter(((t[0][0], t[0][1]), (t[12][0], t[12][1])) for t in trails)
    check("352 input differences, 352 output differences, 6144 pairs",
          len(ins) == 352 and len(outs) == 352 and len(pairs) == 6144,
          f"{len(ins)}/{len(outs)}/{len(pairs)}")
    check("every pair is realized by exactly two characteristics",
          set(pairs.values()) == {2}, f"multiplicities {dict(collections.Counter(pairs.values()))}")
    memo, census = {}, collections.Counter()
    for t in trails:
        v = []
        for r in range(1, 13):
            key = (t[r - 1][0], t[r - 1][1], t[r][0])
            if key not in memo:
                c = exact_count(*key)
                memo[key] = None if (c == 0 or (c & (c - 1))) else 32 - c.bit_length() + 1
            v.append(memo[key])
        census[tuple(v)] += 1
    check("all 12288 have exact weight 58",
          all(None not in v and sum(v) == 58 for v in census), "")
    check("exactly four round-weight vectors, 5376/5376/768/768",
          sorted(census.values(), reverse=True) == [5376, 5376, 768, 768],
          str(sorted(census.values(), reverse=True)))
    hit = [t for t in trails if (t[0][0], t[0][1]) == (DX4, DY4)
           and (t[12][0], t[12][1]) == (A, B)]
    check("the attack's differential is realized by exactly two characteristics, "
          "so its probability is 2 x 2^-58 = 2^-57", len(hit) == 2, str(len(hit)))

print()
print("=" * 78)
print("3.  insertion and the chosen-ciphertext pairing (Step 1)")
print("=" * 78)
check("plaintext-side end is the trail input (forward insertion)",
      (DX4, DY4) == TRAIL[0])
check("ciphertext-side end is the trail output", (A, B) == TRAIL[12])
pair_ok = inv_ok = True
for _ in range(3000):
    X17, Y17, sk16 = rnd(), rnd(), rnd()
    X17s = X17 ^ A
    Y17s = Y17 ^ B ^ (F(X17) ^ F(X17 ^ A))
    X16, Y16 = dec_last(X17, Y17, sk16)
    X16s, Y16s = dec_last(X17s, Y17s, sk16)
    pair_ok &= (X16 ^ X16s == A) and (Y16 ^ Y16s == B)
    X17ss = X17s ^ A
    Y17ss = Y17s ^ B ^ (F(X17s) ^ F(X17s ^ A))
    inv_ok &= (X17ss, Y17ss) == (X17, Y17)
check("pairing gives (dX^16,dY^16) = (a,b) for every sk^16", pair_ok)
check("pairing is a fixed-point-free involution "
      "(2^62 pairs from 2^63 ciphertexts)", inv_ok)

print()
print("=" * 78)
print("4.  key-free filter and the staged conditions")
print("=" * 78)
check("wt(mask of the key-free filter) = 8", bin(FM).count("1") == 8)
check("wt(Stage-1 mask) = 23", bin(S1M).count("1") == 23)
and0 = and1 = M32
or0 = or1 = 0
f_ok = True
for _ in range(20000):
    X4, Y4 = rnd(), rnd()
    sk = [rnd() for _ in range(4)]
    x, y, xs, ys = X4, Y4, X4 ^ DX4, Y4 ^ DY4
    for r in (3, 2, 1, 0):
        x, y = dec(x, y, sk[r])
        xs, ys = dec(xs, ys, sk[r])
    dX0, dY0 = x ^ xs, y ^ ys
    dX1 = dY0 ^ F(x) ^ F(xs)                     # key-independent
    f_ok &= (dX1 & FM) == FV
    and1 &= dX1
    or1 |= dX1
    and0 &= dX0
    or0 |= dX0
check("the key-free filter holds for all 20000 right pairs", f_ok)
c1 = ~(and1 ^ or1) & M32
c0 = ~(and0 ^ or0) & M32
check("its mask is exactly the observed constant bits of dX^1",
      (FM & ~c1) == 0, f"observed {bin(c1).count('1')}, claimed 8")
check("no bit of dX^0 is determined, as the paper states",
      bin(c0).count("1") == 0, f"observed {bin(c0).count('1')}")

# the four incremental conditions must exhaust the 64-bit condition at round 4
stage_ok = True
d3, d4 = collections.Counter(), collections.Counter()
for _ in range(40):
    X4, Y4 = rnd(), rnd()
    sk = [rnd() for _ in range(4)]
    x, y, xs, ys = X4, Y4, X4 ^ DX4, Y4 ^ DY4
    for r in (3, 2, 1, 0):
        x, y = dec(x, y, sk[r])
        xs, ys = dec(xs, ys, sk[r])
    X0, Y0, X0s, Y0s = x, y, xs, ys
    dX0 = X0 ^ X0s
    dX1 = (Y0 ^ Y0s) ^ F(X0) ^ F(X0s)
    stage_ok &= (dX1 & FM) == FV
    X1, X1s = Y0 ^ F(X0) ^ sk[0], Y0s ^ F(X0s) ^ sk[0]
    dX2 = dX0 ^ F(X1) ^ F(X1s)
    stage_ok &= (dX2 & S1M) == S1V
    a, o = M32, 0                                # bits of dX^3 fixed before sk^1
    for _ in range(300):
        g = rnd()
        X2, X2s = X0 ^ F(X1) ^ g, X0s ^ F(X1s) ^ g
        v = dX1 ^ F(X2) ^ F(X2s)
        a &= v
        o |= v
    d3[bin(~(a ^ o) & M32).count("1")] += 1
    a, o = M32, 0                                # bits of dX^4 fixed before sk^2
    for _ in range(300):
        X3 = rnd()
        v = dX2 ^ F(X3) ^ F(X3 ^ DY4)
        a &= v
        o |= v
    d4[bin(~(a ^ o) & M32).count("1")] += 1
check("right pairs satisfy the key-free filter and the Stage-1 condition", stage_ok)
check("the Stage-3 test is a 9-bit condition (23 bits of dX^4 are pre-fixed)",
      set(d4) == {23}, f"observed {dict(sorted(d4.items()))}")
print(f"  bits of dX^3 pre-fixed before guessing sk^1: {dict(sorted(d3.items()))}")
print(f"  (varies per pair; the average filtering factor is 2^-24, since")
print(f"   8 + 23 + 24 + 9 = 64 exhausts the condition at round 4)")
check("the four conditions sum to the 64-bit condition", 8 + 23 + 24 + 9 == 64)

print()
print("=" * 78)
print("5.  guessed-key masks (Appendix A)")
print("=" * 78)
w = lambda x: bin(x).count("1")
check("wt(M_K0) = 32, wt(M_K1) = 28, wt(M_K2) = 9",
      w(MK0) == 32 and w(MK1) == 28 and w(MK2) == 9,
      f"{w(MK0)}/{w(MK1)}/{w(MK2)}")
check("M_K1 = M_K1^(1) OR M_K1^(2), weights 22 and 20",
      MK1 == (MK1a | MK1b) and w(MK1a) == 22 and w(MK1b) == 20)
check("k = 32 + 28 + 9 = 69", w(MK0) + w(MK1) + w(MK2) == 69)
check("remaining master-key bits = 128 - 69 = 59", 128 - 69 == 59)
det = pre = True
for _ in range(3000):
    X3 = rnd()
    base = F(X3) ^ F(X3 ^ DY4)                   # dX^3 = dY^4 = 008C0000
    for _ in range(4):
        X3b = (X3 & MK2) | (rnd() & ~MK2 & M32)
        det &= (F(X3b) ^ F(X3b ^ DY4)) == base
    X2 = rnd()
    b2 = F(X2) & MK2
    for _ in range(4):
        X2b = (X2 & MK1b) | (rnd() & ~MK1b & M32)
        pre &= (F(X2b) & MK2) == b2
check("D_{dX^3} f(X^3) depends only on the 9 M_K2 coordinates", det)
check("f(X^2) restricted to M_K2 depends only on the M_K1^(2) coordinates", pre)
comp = True
for _ in range(400):
    X4, Y4 = rnd(), rnd()
    sk = [rnd() for _ in range(4)]
    x, y, xs, ys = X4, Y4, X4 ^ DX4, Y4 ^ DY4
    for r in (3, 2, 1, 0):
        x, y = dec(x, y, sk[r])
        xs, ys = dec(xs, ys, sk[r])
    X0, Y0, X0s, Y0s = x, y, xs, ys

    def test(g0, g1, g2):
        X1, X1s = Y0 ^ F(X0) ^ g0, Y0s ^ F(X0s) ^ g0
        X2, X2s = X0 ^ F(X1) ^ g1, X0s ^ F(X1s) ^ g1
        X3, X3s = X1 ^ F(X2) ^ g2, X1s ^ F(X2s) ^ g2
        return (X2 ^ X2s, X1 ^ X1s ^ F(X2) ^ F(X2s), X2 ^ X2s ^ F(X3) ^ F(X3s))

    ref = test(sk[0], sk[1], sk[2])
    comp &= ref[1] == DY4 and ref[2] == DX4      # dX^3 and dX^4 as required
    for _ in range(6):
        comp &= test((sk[0] & MK0) | (rnd() & ~MK0 & M32),
                     (sk[1] & MK1) | (rnd() & ~MK1 & M32),
                     (sk[2] & MK2) | (rnd() & ~MK2 & M32)) == ref
check("right pairs are counted by the true key class, and the test is invariant "
      "outside M_K0/M_K1/M_K2", comp)

print()
print("=" * 78)
print("6.  complexity arithmetic and success probability")
print("=" * 78)
NP, WT, k, rounds = 62, 57, 69, 17               # pairs, differential weight, key bits
t = [2 ** NP, 2 ** 32 * 2 ** 54, 2 ** (32 + 22) * 2 ** 31, 2 ** (32 + 22 + 15) * 2 ** 7]
for n, v in zip(["Step 2 (key-free filter)", "Stage 1", "Stage 2", "Stage 3"], t):
    print(f"  {n:<26} = 2^{math.log2(v):>5.0f}")
core = (2 / rounds) * sum(t)
T = math.log2(core + 2 ** (128 - k))
M = math.log2(2 ** 54 * 4 * 32 / 8)
print(f"  (2/17) x sum + 2^59 = 2^{T:.2f}   (paper 2^83.5)")
print(f"  M = 2^{M:.0f} bytes                (paper 2^58)")
check("time matches the paper", abs(T - 83.5) < 0.1, f"2^{T:.2f}")
check("memory matches the paper", abs(M - 58) < 0.1, f"2^{M:.2f}")
check("data 2^63 < full codebook 2^64", NP + 1 < 64)

try:
    from decimal import Decimal as D, getcontext
    getcontext().prec = 90
    pmf = lambda i, lam: (-lam).exp() * lam ** i / D(math.factorial(i))
    lam_r, lam_w, m = D(2) ** (NP - WT), D(2) ** (NP - 64), D(2) ** k - 1
    succ = 1 - sum(pmf(i, lam_r) for i in range(16))
    kept = m * (1 - sum(pmf(i, lam_w) for i in range(16)))
    print(f"  E[c_right] = {float(lam_r):g}, E[c_wrong] = {float(lam_w):g}, "
          f"classes 2^{k}")
    print(f"  tau = 16: P(success) = {float(succ):.5f}, "
          f"E[# wrong classes retained] = 2^{math.log2(float(kept)):.1f}")
    check("success probability matches the paper (0.99934)",
          abs(float(succ) - 0.99934) < 5e-5, f"{float(succ):.5f}")
    check("no wrong class is expected to survive tau = 16", float(kept) < 0.01,
          f"{float(kept):.4f}")
except ImportError:
    print("  [skipped] decimal not available")

print()
print("=" * 78)
print("7.  minimum differential weights B[r] (Table 5)")
print("=" * 78)
bpath = os.path.join(HERE, "bounds_sand64.json")
if not os.path.exists(bpath):
    print(f"  [skipped] {bpath} not found")
else:
    Bv = {int(a): b for a, b in json.load(open(bpath)).items()}
    paper = {1: 0, 2: 2, 3: 4, 4: 8, 5: 12, 6: 18, 7: 24, 8: 32,
             9: 39, 10: 45, 11: 52, 12: 58, 13: 64, 14: 68, 15: 72}
    check("B[r] agrees with Table 5 for every r <= 15",
          all(Bv.get(r) == v for r, v in paper.items()), str(Bv))
    a_des = {2: 1, 3: 2, 4: 4, 5: 6, 6: 9, 7: 12, 8: 16,
             9: 18, 10: 20, 11: 22, 12: 24, 13: 26, 14: 29, 15: 32}
    eq = [r for r in a_des if Bv[r] == 2 * a_des[r]]
    gt = [r for r in a_des if Bv[r] > 2 * a_des[r]]
    check("B[r] = 2a[r] exactly for r <= 8 and B[r] > 2a[r] from r = 9 on",
          eq == list(range(2, 9)) and gt == list(range(9, 16)),
          f"equal {eq}, strict {gt}")
    check("B[13] = 64, so the 2^-64 level is reached at 13 rounds "
          "(designers: 15)", Bv[13] == 64)

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
print("=" * 78)
sys.exit(0 if ok_all else 1)
