#!/usr/bin/env python3
"""Random-difference heuristic for the 16-round chosen-ciphertext attack.

Chosen-ciphertext pairs with (dx15,dy15) = (00182C28,00404000) are decrypted
through 16 rounds (independent random round keys).  We measure the pass rate
of Eq. (filter0) (8-bit mask, model 2^-8) and, separately, the number of
pairs passing Eq. (filter1) (23-bit mask, model 2^-23).
"""
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

R = 16
sk = [int(v) for v in rnd(R)]
A, B = 0x00182C28, 0x00404000
M0F, V0 = 0x81311081, 0x00110000
M1F, V1 = 0xB9FF37B9, 0x001E0000
tot0 = tot1 = totb = n = 0
for _ in range(16):
    N = 1 << 22
    c0, c1 = rnd(N), rnd(N)
    d0 = c0 ^ U32(A)
    d1 = c1 ^ U32(B) ^ f(c0) ^ f(d0)
    p = dec(c0, c1, sk); q = dec(d0, d1, sk)
    dx0 = p[0] ^ q[0]
    dx1 = (p[1] ^ q[1]) ^ f(p[0]) ^ f(q[0])
    f0 = (dx0 & U32(M0F)) == U32(V0)
    f1 = (dx1 & U32(M1F)) == U32(V1)
    tot0 += int(f0.sum()); tot1 += int(f1.sum()); totb += int((f0 & f1).sum()); n += N
print(f"16 rounds, {n} = 2^{np.log2(n):.0f} chosen-ciphertext pairs")
print(f"  filter0 rate 2^{np.log2(tot0/n):.3f} (model 2^-8), count {tot0}")
print(f"  filter1 count {tot1} (model {n*2.0**-23:.1f}), both {totb} (model {n*2.0**-31:.4f})")
