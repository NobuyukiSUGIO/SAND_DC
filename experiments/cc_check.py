#!/usr/bin/env python3
"""Experimental check of the chosen-ciphertext attack of the paper (Sect. V-C).

16-round attack: rounds 0-2 key recovery (sk^0, sk^1), rounds 3-14 the
characteristic, round 15 swap-free.
Toy version (8 rounds): the same three key-recovery rounds, a 4-round
characteristic in the forward (encryption) direction
    (dx3,dy3) = (00900000,008C0000) -> (dx7,dy7) = (00900000,00010000)
(weight 4+2+0+2 = 8), and the swap-free last round.  The output difference
(dx7,dy7) of the characteristic is input from the ciphertext side: ciphertext
pairs are chosen so that (dx7,dy7) holds independently of sk^7, and the
plaintexts are obtained from the decryption oracle.

Checks:
 [1] chosen-ciphertext pairs give (dx7,dy7) independently of sk^7
 [2] right pairs (reaching (dx3,dy3)) pass Eqs. (filter0) and (filter1);
     pass rates of the filters
 [3] the right class counts every right pair; the test depends only on
     sk0 & DFDFDFDF and sk1 & 8C008C8C
 [4] near-key classes: fraction phi of right pairs they count; expected
     number of wrong classes with counter >= tau = 16 for lambda = 32, and
     the probability that the right class is not ranked first
"""
import itertools, math, random
import numpy as np
# --- SAND-64 implementation (row-wise words: x = x{3}||x{2}||x{1}||x{0},
# --- x{i} = bits 8i+7..8i), self-contained -------------------------------
rng = np.random.default_rng(2026)
U32 = np.uint32
P8 = [7, 4, 1, 6, 3, 0, 5, 2]


def word(x, i):
    return (x >> U32(8 * i)) & U32(0xFF)


def pack(w3, w2, w1, w0):
    return (w3 << U32(24)) | (w2 << U32(16)) | (w1 << U32(8)) | w0


def G0(x):
    x3, x2, x1, x0 = word(x, 3), word(x, 2), word(x, 1), word(x, 0)
    y0 = (x3 & x2) ^ x0
    y3 = (y0 & x1) ^ x3
    return pack(y3, x2, x1, y0)


def G1(x):
    x3, x2, x1, x0 = word(x, 3), word(x, 2), word(x, 1), word(x, 0)
    y2 = (x3 & x1) ^ x2
    y1 = (y2 & x0) ^ x1
    return pack(x3, y2, y1, x0)


def rotw(x, t):
    out = U32(0) * x
    for i in range(4):
        w = word(x, i)
        w = ((w << U32(t)) | (w >> U32(8 - t))) & U32(0xFF)
        out |= w << U32(8 * i)
    return out


def P32(x):
    out = U32(0) * x
    for i in range(4):
        for j in range(8):
            out |= ((x >> U32(8 * i + j)) & U32(1)) << U32(8 * i + P8[j])
    return out


def f(x):
    return P32(G0(x) ^ G1(rotw(x, 1)))


def enc(xl, yr, sk):
    x, y = xl.copy(), yr.copy()
    R = len(sk)
    for r in range(R - 1):
        x, y = y ^ f(x) ^ U32(sk[r]), x
    return x, y ^ f(x) ^ U32(sk[R - 1])


def dec(xl, yr, sk):
    # Feistel without final swap: decryption = encryption with reversed keys
    return enc(xl, yr, sk[::-1])


def rnd(n):
    return rng.integers(0, 2**32, size=n, dtype=np.uint64).astype(U32)

R = 8
sk = [int(v) for v in rnd(R)]
A, B = 0x00900000, 0x00010000               # (dx7, dy7) in the toy
TX3, TY3 = 0x00900000, 0x008C0000
M0F, V0 = 0x81311081, 0x00110000
M1F, V1 = 0xB9FF37B9, 0x001E0000
MK0, MK1 = 0xDFDFDFDF, 0x8C008C8C

# [0] the 4-round characteristic holds in the forward direction ---------------
TRAIL4 = [(0x00900000, 0x008C0000), (0x00010000, 0x00900000),
          (0x00000000, 0x00010000), (0x00010000, 0x00000000),
          (0x00900000, 0x00010000)]          # (dx^r, dy^r), r = 3..7
xs = rnd(1 << 22)
ws = []
for r in range(4):
    (dx, dy), (dxn, dyn) = TRAIL4[r], TRAIL4[r + 1]
    assert dyn == dx
    ws.append(-np.log2(np.mean((f(xs) ^ f(xs ^ U32(dx))) == U32(dy ^ dxn))))
print("[0] forward round weights r=3..6:", [round(w, 2) for w in ws], "total", round(sum(ws), 2))

# [1] chosen ciphertext pairs ------------------------------------------------
Np = 1 << 20
c0, c1 = rnd(Np), rnd(Np)
d0 = c0 ^ U32(A)
d1 = c1 ^ U32(B) ^ f(c0) ^ f(d0)
p = dec(c0, c1, sk)
q = dec(d0, d1, sk)


def fwd(x, y, ks):
    for k in ks:
        x, y = y ^ f(x) ^ U32(k), x
    return x, y


s7p, s7q = fwd(*p, sk[:R - 1]), fwd(*q, sk[:R - 1])
ok7 = np.all((s7p[0] ^ s7q[0]) == U32(A)) and np.all((s7p[1] ^ s7q[1]) == U32(B))
print(f"[1] {Np} chosen-ciphertext pairs, (dx7,dy7) = (00900000,00010000) for all: {bool(ok7)}")

# [2] right pairs and filters ---------------------------------------------------
s3p, s3q = fwd(*p, sk[:3]), fwd(*q, sk[:3])
right = ((s3p[0] ^ s3q[0]) == U32(TX3)) & ((s3p[1] ^ s3q[1]) == U32(TY3))
dx0 = p[0] ^ q[0]
dx1 = (p[1] ^ q[1]) ^ f(p[0]) ^ f(q[0])            # key independent
f0 = (dx0 & U32(M0F)) == U32(V0)
f1 = f0 & ((dx1 & U32(M1F)) == U32(V1))
print(f"[2] right pairs {right.sum()} (rate 2^{np.log2(right.mean()):.2f}, characteristic 2^-8)")
print(f"    filter0 pass {f0.sum()}, filter0+1 pass {f1.sum()}; all right pairs pass: {bool(np.all(f1[right]))}")
nr = ~right
print(f"    among non-right pairs: filter0 rate 2^{np.log2(f0[nr].mean()):.2f} (model 2^-8)")

# [3] key test ---------------------------------------------------------------
idx = np.nonzero(right)[0]
P = (p[0][idx], p[1][idx], q[0][idx], q[1][idx])


def kr_test(px, py, qx, qy, k0, k1):
    x1p = py ^ f(px) ^ U32(k0); x1q = qy ^ f(qx) ^ U32(k0)
    x2p = px ^ f(x1p) ^ U32(k1); x2q = qx ^ f(x1q) ^ U32(k1)
    return ((x2p ^ x2q) == U32(TY3)) & (((x1p ^ x1q) ^ f(x2p) ^ f(x2q)) == U32(TX3))


base = kr_test(*P, sk[0], sk[1])
print(f"[3] right class counts {base.sum()} of {len(idx)} right pairs")
same = True
for _ in range(200):
    k0 = (sk[0] & MK0) | (int(rnd(1)[0]) & ~MK0 & 0xFFFFFFFF)
    k1 = (sk[1] & MK1) | (int(rnd(1)[0]) & ~MK1 & 0xFFFFFFFF)
    same &= bool(np.all(kr_test(*P, k0, k1) == base))
print("    test invariant under key bits outside DFDFDFDF / 8C008C8C (200 trials):", same)

# [4] near-key classes ----------------------------------------------------------
bits = [(0, i) for i in range(32) if MK0 >> i & 1] + [(1, i) for i in range(32) if MK1 >> i & 1]


def phi(flips):
    k = [sk[0], sk[1]]
    for t, i in flips:
        k[t] ^= 1 << i
    return kr_test(*P, k[0], k[1]).mean()


def pois_pmf(k, lam): return math.exp(-lam) * lam**k / math.factorial(k)
def pois_sf(k, lam): return 1 - sum(pois_pmf(i, lam) for i in range(k + 1))
def binom_pmf(k, n, x): return math.comb(n, k) * x**k * (1 - x)**(n - k)


LAM, TAU, NOISE = 32.0, 16, 0.25
random.seed(1)
tot_ret = tot_rank = 0.0
for d in range(1, 9):
    n_all = math.comb(len(bits), d)
    if n_all <= 8000:
        combos, scale = list(itertools.combinations(bits, d)), 1.0
    else:
        combos, scale = [random.sample(bits, d) for _ in range(8000)], n_all / 8000
    ph = [x for x in (phi(list(c)) for c in combos) if x > 0]
    ret = scale * sum(pois_sf(TAU - 1, LAM * x + NOISE) for x in ph)
    rank = 0.0
    for x in ph:
        sacc = 0.0
        for c in range(1, 80):
            qv = sum(binom_pmf(k, c, x) * (pois_sf(c - k - 1, NOISE) if c - k - 1 >= 0 else 1.0)
                     for k in range(c + 1))
            sacc += pois_pmf(c, LAM) * qv
        rank += scale * sacc
    tot_ret += ret; tot_rank += rank
    print(f"[4] d={d}: max phi {max(ph) if ph else 0:.3f}, E[#retained]={ret:.3f}, P(outrank)<={rank:.2e}")
print(f"[4] total: E[#wrong retained] ~ {tot_ret:.1f}, P(right not ranked first) <= {tot_rank:.2e}")
print(f"    Pr[Pois(32) >= 16] = {pois_sf(15, 32.0):.5f}")
