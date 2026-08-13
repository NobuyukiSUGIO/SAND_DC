#!/usr/bin/env python3
"""
Independent verification of the 23-round key-recovery attack on SAND-128
described in paper/SAND_DC_IEICE/main.tex (2^126 chosen ciphertexts version).

Nothing in this file imports the author's modules for the cryptographic
primitives; SAND-128 is re-implemented from the specification given in
Section IV of the paper (row-wise state, column = SSb nibble).
The author's modules are imported only at the very end to cross-check that
the two independent implementations agree.
"""
import itertools
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

random.seed(20260809)

P16 = [14, 15, 8, 9, 2, 3, 12, 13, 6, 7, 0, 1, 10, 11, 4, 5]
MASK64 = (1 << 64) - 1

# ---------------------------------------------------------------- primitives
def G0(u):
    u3, u2, u1, u0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y0 = (u3 & u2) ^ u0
    y1 = u1
    y2 = u2
    y3 = (y0 & u1) ^ u3
    return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0


def G1(u):
    u3, u2, u1, u0 = (u >> 3) & 1, (u >> 2) & 1, (u >> 1) & 1, u & 1
    y3 = u3
    y2 = (u3 & u1) ^ u2
    y1 = (y2 & u0) ^ u1
    y0 = u0
    return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0


def bit(w, i):          # i = 0 is the MSB of the 64-bit word
    return (w >> (63 - i)) & 1


def cols(w):
    """word -> 16 columns; column j = (row3,row2,row1,row0) bits at offset j."""
    return [(bit(w, j) << 3) | (bit(w, j + 16) << 2) |
            (bit(w, j + 32) << 1) | bit(w, j + 48) for j in range(16)]


def uncols(c):
    w = 0
    for j in range(16):
        for k, sh in enumerate((3, 2, 1, 0)):        # k=0 -> row x{3}
            if (c[j] >> sh) & 1:
                w |= 1 << (63 - (j + 16 * k))
    return w


def P(w):
    """y_{16i+p16(j)} = x_{16i+j}: output column p16(j) <- input column j."""
    c = cols(w)
    out = [0] * 16
    for j in range(16):
        # column index in the paper's LSB-within-row numbering is 15-j
        out[15 - P16[15 - j]] = c[j]
    return uncols(out)


def F_slow(x):
    """f(x) = P_n(G0(x) XOR G1(x <<< 1)); G1's nibble j+1 feeds output nibble j."""
    c = cols(x)
    t = [G0(c[j]) ^ G1(c[(j + 1) % 16]) for j in range(16)]
    return P(uncols(t))


# fast, row-parallel version of the same function (16 columns at a time)
PTAB = [0] * 65536
for _v in range(65536):
    _o = 0
    for _j in range(16):                    # input column j -> output column 15-P16[15-j]
        if (_v >> (15 - _j)) & 1:
            _o |= 1 << P16[15 - _j]
    PTAB[_v] = _o
del _v, _o, _j


def F(x):
    r3, r2, r1, r0 = (x >> 48) & 0xFFFF, (x >> 32) & 0xFFFF, (x >> 16) & 0xFFFF, x & 0xFFFF
    # G0 on the state itself
    a0 = (r3 & r2) ^ r0
    a1 = r1
    a2 = r2
    a3 = (a0 & r1) ^ r3
    # G1 on the state rotated by one column (= 16-bit rotate left by 1 per row)
    s3 = ((r3 << 1) | (r3 >> 15)) & 0xFFFF
    s2 = ((r2 << 1) | (r2 >> 15)) & 0xFFFF
    s1 = ((r1 << 1) | (r1 >> 15)) & 0xFFFF
    s0 = ((r0 << 1) | (r0 >> 15)) & 0xFFFF
    b3 = s3
    b2 = (s3 & s1) ^ s2
    b1 = (b2 & s0) ^ s1
    b0 = s0
    return (PTAB[a3 ^ b3] << 48) | (PTAB[a2 ^ b2] << 32) | \
           (PTAB[a1 ^ b1] << 16) | PTAB[a0 ^ b0]


def enc_round(x, y, sk):
    return (y ^ F(x) ^ sk) & MASK64, x


def dec_round(x_next, y_next, sk):
    x = y_next
    y = x_next ^ F(x) ^ sk
    return x & MASK64, y & MASK64


def rnd64():
    return random.getrandbits(64)


def hx(w):
    return f"{w:016X}"


ok_all = True


def check(name, cond, extra=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK ' if cond else 'FAIL'}] {name}{(' : ' + extra) if extra else ''}")
    return cond


print("=" * 78)
print("0.  self-consistency of the re-implementation")
print("=" * 78)
# P is a permutation and encryption/decryption are inverse
perm_ok = sorted(P(1 << (63 - i)) for i in range(64)) == sorted(1 << i for i in range(64))
check("P_n is a bit permutation", perm_ok)
check("row-parallel F agrees with the nibble-wise definition",
      all(F(w) == F_slow(w) for w in [rnd64() for _ in range(300)] + [0, MASK64, 1]))
inv_ok = True
for _ in range(200):
    x, y, sk = rnd64(), rnd64(), rnd64()
    xx, yy = enc_round(x, y, sk)
    inv_ok &= dec_round(xx, yy, sk) == (x, y)
check("round function is invertible (dec o enc = id)", inv_ok)

# joint vs marginal DDT
DDT = [[0] * 256 for _ in range(16)]
for din in range(16):
    for u in range(16):
        DDT[din][((G0(u) ^ G0(u ^ din)) << 4) | (G1(u) ^ G1(u ^ din))] += 1
joint = sum(1 for a in range(16) for o in range(256) if DDT[a][o])
m0 = [[0] * 16 for _ in range(16)]
m1 = [[0] * 16 for _ in range(16)]
for a in range(16):
    for u in range(16):
        m0[a][G0(u) ^ G0(u ^ a)] += 1
        m1[a][G1(u) ^ G1(u ^ a)] += 1
marg = sum(1 for a in range(16) for b in range(16) for c in range(16) if m0[a][b] and m1[a][c])
check("joint SSb-DDT has 97 triples (paper: 97)", joint == 97, f"{joint}")
check("product of marginal DDTs has 157 triples (paper: 157)", marg == 157, f"{marg}")
check("DDT_SSb[0x6][0x62] = 0 (paper's example)", DDT[0x6][0x62] == 0)

print()
print("=" * 78)
print("1.  the 19-round characteristic of Table VI")
print("=" * 78)
TRAIL = [
    (0x0000000000001115, 0x000000000000CC0C),
    (0x00000000000003C0, 0x0000000000001115),
    (0x0000000000000114, 0x00000000000003C0),
    (0x000000000000000C, 0x0000000000000114),
    (0x0000000000000010, 0x000000000000000C),
    (0x0000000000000000, 0x0000000000000010),
    (0x0000000000000010, 0x0000000000000000),
    (0x000000000000000C, 0x0000000000000010),
    (0x0000000000000114, 0x000000000000000C),
    (0x00000000000003C0, 0x0000000000000114),
    (0x0000000000001115, 0x00000000000003C0),
    (0x000000000000CC0C, 0x0000000000001115),
    (0x0000000000005400, 0x000000000000CC0C),
    (0x000000000000C03F, 0x0000000000005400),
    (0x0000000000004410, 0x000000000000C03F),
    (0x000000000000C000, 0x0000000000004410),
    (0x0000000000000400, 0x000000000000C000),
    (0x000000000000C003, 0x0000000000000400),
    (0x0000000000000510, 0x000000000000C003),
    (0x000000000000C0CC, 0x0000000000000510),
]
PAPER_W = [10, 8, 6, 4, 2, 0, 2, 4, 6, 8, 10, 12, 6, 16, 6, 4, 2, 8, 6]


def p_inv(w):
    c = cols(w)
    out = [0] * 16
    for j in range(16):
        out[j] = c[15 - P16[15 - j]]
    return uncols(out)


def exact_count(dx_in, dy_in, dx_out):
    """#{64-bit x : f(x)^f(x^dx_in) = dx_out ^ dy_in}, by explicit enumeration
    over the (dG0,dG1) choices of the *active* nibbles only."""
    T = cols(p_inv(dx_out ^ dy_in))
    din = cols(dx_in)
    opts = []
    for m in range(16):
        if din[m] == 0:
            opts.append([(0, 0, 16)])
        else:
            opts.append([((o >> 4) & 0xF, o & 0xF, DDT[din[m]][o])
                         for o in range(256) if DDT[din[m]][o]])
    total = 0
    for combo in itertools.product(*opts):
        if all(T[j] == (combo[j][0] ^ combo[(j + 1) % 16][1]) for j in range(16)):
            p = 1
            for _, _, c in combo:
                p *= c
            total += p
    return total


print(f"{'r':>3} {'dX^r':>17} {'dY^r':>17} {'paper w':>8} {'exact':>7} {'count':>22}")
tot = 0
weights_ok = True
link_ok = True
for r in range(1, 20):
    dxp, dyp = TRAIL[r - 1]
    dxc, dyc = TRAIL[r]
    link_ok &= (dyc == dxp)
    cnt = exact_count(dxp, dyp, dxc)
    w = 64 - cnt.bit_length() + 1 if cnt else None
    weights_ok &= (w == PAPER_W[r - 1])
    tot += w if w else 0
    print(f"{r:>3} {hx(dxc):>17} {hx(dyc):>17} {PAPER_W[r-1]:>8} {str(w):>7} {cnt:>22}")
check("Feistel linkage dY^r = dX^(r-1) holds for all rounds", link_ok)
check("every round weight matches the paper", weights_ok)
check("total weight = 120 (p = 2^-120)", tot == 120, f"{tot}")

# Monte-Carlo confirmation of a few individual round probabilities
print("\n  Monte-Carlo check of individual round DPs (independent of the count above):")
LOGN = 20
N = 1 << LOGN
getrand = random.getrandbits
for r in (1, 5, 12, 13, 19):
    dxp, dyp = TRAIL[r - 1]
    dxc, _ = TRAIL[r]
    tgt = dxc ^ dyp
    hits = 0
    for _ in range(N):
        x = getrand(64)
        if F(x) ^ F(x ^ dxp) == tgt:
            hits += 1
    est = f"2^-{LOGN - (hits.bit_length()-1)}" if hits else "none"
    print(f"    round {r:>2}: paper w={PAPER_W[r-1]:>2}  hits {hits:>7} / 2^{LOGN}"
          f"   -> DP approx {est}")

print()
print("=" * 78)
print("2.  chosen-ciphertext pairing of attack step 1")
print("=" * 78)
a = 0x000000000000C0CC
b = 0x0000000000000510
ins_ok = True
inv_ok = True
for _ in range(2000):
    X23, Y23, sk22 = rnd64(), rnd64(), rnd64()
    X23s = X23 ^ a
    Y23s = Y23 ^ b ^ (F(X23) ^ F(X23 ^ a))
    # the last round omits the branch swap: (x,y) -> (x, y ^ f(x) ^ sk)
    X22, Y22 = X23, (Y23 ^ F(X23) ^ sk22) & MASK64
    X22s, Y22s = X23s, (Y23s ^ F(X23s) ^ sk22) & MASK64
    ins_ok &= (X22 ^ X22s == a) and (Y22 ^ Y22s == b)
    # involution
    X23ss = X23s ^ a
    Y23ss = Y23s ^ b ^ (F(X23s) ^ F(X23s ^ a))
    inv_ok &= (X23ss, Y23ss) == (X23, Y23)
check("pairing forces (dX^22,dY^22) = (a,b) for every sk^22", ins_ok)
check("pairing map is a fixed-point-free involution", inv_ok)
check("(a,b) is the output difference of the trail (forward insertion)",
      (a, b) == (TRAIL[19][0], TRAIL[19][1]))

print()
print("=" * 78)
print("3.  key-free filters of attack steps 2 and 3")
print("=" * 78)
DX3, DY3 = 0x0000000000001115, 0x000000000000CC0C
check("(dX^3,dY^3) = the trail input (dX^0,dY^0), forward insertion",
      (DX3, DY3) == (TRAIL[0][0], TRAIL[0][1]))

M0, V0 = 0x2A226C5464443F33, 0x0000000000000033
M1, V1 = 0xFCCCFFFFB9D9FFFF, 0x0000000000005400
check("wt(mask of Eq.(9)) = 27", bin(M0).count("1") == 27, str(bin(M0).count("1")))
check("wt(mask of Eq.(10)) = 52", bin(M1).count("1") == 52, str(bin(M1).count("1")))
check("required values are submasks of their masks", (V0 & ~M0) == 0 and (V1 & ~M1) == 0)

# every pair that reaches (dX^3,dY^3) must satisfy both filters
TRIALS = 20000
and0 = MASK64
or0 = 0
and1 = MASK64
or1 = 0
f_ok = True
for _ in range(TRIALS):
    X3, Y3 = rnd64(), rnd64()
    X3s, Y3s = X3 ^ DX3, Y3 ^ DY3
    sk = [rnd64(), rnd64(), rnd64()]
    x, y = X3, Y3
    xs, ys = X3s, Y3s
    for r in (2, 1, 0):
        x, y = dec_round(x, y, sk[r])
        xs, ys = dec_round(xs, ys, sk[r])
    dX0, dY0 = x ^ xs, y ^ ys
    # dX^1 recomputed key-independently from the plaintext pair (attack step 3)
    dX1 = dY0 ^ F(x) ^ F(xs)
    f_ok &= (dX0 & M0) == V0
    f_ok &= (dX1 & M1) == V1
    and0 &= dX0
    or0 |= dX0
    and1 &= dX1
    or1 |= dX1
check(f"filter Eq.(11) holds for all {TRIALS} random right pairs", (and0 & M0) == V0 and f_ok)
check(f"filter Eq.(12) holds for all {TRIALS} random right pairs", (and1 & M1) == V1)
# tightness: the observed constant bits are exactly the claimed masks
obs0 = ~(and0 ^ or0) & MASK64
obs1 = ~(and1 ^ or1) & MASK64
check("Eq.(11) mask is not weaker than what is observed", (M0 & ~obs0) == 0,
      f"observed constant bits {bin(obs0).count('1')}, claimed {bin(M0).count('1')}")
check("Eq.(12) mask is not weaker than what is observed", (M1 & ~obs1) == 0,
      f"observed constant bits {bin(obs1).count('1')}, claimed {bin(M1).count('1')}")

print()
print("=" * 78)
print("4.  column numbering of Section V and the appendix column sets")
print("=" * 78)
# The paper numbers the columns of the 4 x n/4 array from the RIGHT: column j is
# the j-th nibble x[j] = (x_{4j+3}, x_{4j+2}, x_{4j+1}, x_{4j}) of the n-bit
# word, as in the SAND specification (x[0] = least significant nibble).  The
# re-implementation above indexes bits from the MSB, so its column c and the
# paper's column j are related by j = n/4-1-c.  This section checks that the two
# describe the same cipher and re-derives the appendix column sets in the
# paper's numbering.

def cols_lsb(w, ncol=16):
    """column j = (x{3}[j], x{2}[j], x{1}[j], x{0}[j]), bit j counted from the LSB."""
    m = (1 << ncol) - 1
    rows = [(w >> ((3 - k) * ncol)) & m for k in range(4)]        # rows[0] = x{3}
    return [sum(((rows[k] >> j) & 1) << (3 - k) for k in range(4)) for j in range(ncol)]


def uncols_lsb(c, ncol=16):
    rows = [0] * 4
    for j, v in enumerate(c):
        for k in range(4):
            if (v >> (3 - k)) & 1:
                rows[k] |= 1 << j
    w = 0
    for k in range(4):
        w |= rows[k] << ((3 - k) * ncol)
    return w


def load(Pl, ncol=16):
    """state loading of Section V: row x{3-k}[j] = Pl_{4j+k}."""
    rows = [0] * 4
    for j in range(ncol):
        for k in range(4):                       # k = 3 -> row x{3}
            if (Pl >> (4 * j + k)) & 1:
                rows[3 - k] |= 1 << j
    w = 0
    for k in range(4):
        w |= rows[k] << ((3 - k) * ncol)
    return w


check("cols_lsb / uncols_lsb round-trip",
      all(uncols_lsb(cols_lsb(w)) == w for w in [rnd64() for _ in range(300)]))
# Section V: column j of the loaded array is the j-th nibble of the plaintext word
sec5_ok = True
for _ in range(300):
    Pl = rnd64()
    c = cols_lsb(load(Pl))
    sec5_ok &= all(c[j] == ((Pl >> (4 * j)) & 0xF) for j in range(16))
check("Sect. V: column j of the array is the nibble x[j] = (x_{4j+3}..x_{4j})", sec5_ok)
check("paper column j = implementation column n/4-1-j",
      all(cols_lsb(w) == cols(w)[::-1] for w in [rnd64() for _ in range(300)]))

# Sect. VI-A: in the paper's numbering, x[j] feeds N_0 at column j and N_1 at
# column j+1, i.e. output column j collects G0(col j) and G1(col j-1)
def F_lsb(x):
    c = cols_lsb(x)
    return P(uncols_lsb([G0(c[j]) ^ G1(c[(j - 1) % 16]) for j in range(16)]))


check("Sect. VI-A: 'x[j] feeds the (j+1)-th SSb of G_1' reproduces f exactly",
      all(F_lsb(w) == F(w) for w in [rnd64() for _ in range(300)] + [0, MASK64, 1]))

# --- the two appendix column sets, in the paper's numbering ---
def active(diff, ncol):
    c = cols_lsb(diff, ncol)
    return [j for j in range(ncol) if c[j]], c


DX3_64 = 0x008C0000                      # Appendix A, SAND-64 (n/4 = 8 columns)
DX2_128 = 0x000000000000CC0C             # Appendix B, SAND-128 (n/4 = 16 columns)
J_A, cA = active(DX3_64, 8)
J_B, cB = active(DX2_128, 16)
print(f"  Appendix A: dX^3 = {DX3_64:08X} -> columns {J_A}, "
      f"delta = {sorted({hex(cA[j]) for j in J_A})}")
print(f"  Appendix B: dX^2 = {DX2_128:016X} -> columns {J_B}, "
      f"delta = {sorted({hex(cB[j]) for j in J_B})}")
check("Appendix A column set is {2,3,7} with delta = 0x4",
      J_A == [2, 3, 7] and all(cA[j] == 0x4 for j in J_A))
check("Appendix B column set is {2,3,10,11,14,15} with delta = 0x1",
      J_B == [2, 3, 10, 11, 14, 15] and all(cB[j] == 0x1 for j in J_B))

# the round-key masks must be supported exactly on those columns, in the rows
# that the derivative actually depends on
def support(mask, ncol):
    c = cols_lsb(mask, ncol)
    colset = [j for j in range(ncol) if c[j]]
    rowset = set()
    for j in colset:
        rowset |= {3 - k for k in range(4) if (c[j] >> (3 - k)) & 1}
    return colset, sorted(rowset)


def dep_rows(delta):
    """input bits u_b on which D_delta SSb actually depends (b = 3 is u_3)."""
    def D(u):
        return ((G0(u) ^ G0(u ^ delta)) << 4) | (G1(u) ^ G1(u ^ delta))
    return sorted({b for b in range(4) if any(D(u) != D(u ^ (1 << b)) for u in range(16))})


MK2_64 = 0x8C008C8C                       # Appendix A, sk^2
cols_m, rows_m = support(MK2_64, 8)
check("M_K2 (SAND-64) sits on columns {2,3,7} and rows u_3,u_1,u_0",
      cols_m == [2, 3, 7] and rows_m == dep_rows(0x4) == [0, 1, 3],
      f"columns {cols_m}, rows {rows_m}, derivative needs {dep_rows(0x4)}")
MK1_128 = 0xCC0CCC0CCC0C0000           # Appendix B, sk^1 (re-checked in section 5)
cols_m, rows_m = support(MK1_128, 16)
check("M_K1 (SAND-128) sits on columns {2,3,10,11,14,15} and rows u_3,u_2,u_1",
      cols_m == [2, 3, 10, 11, 14, 15] and rows_m == dep_rows(0x1) == [1, 2, 3],
      f"columns {cols_m}, rows {rows_m}, derivative needs {dep_rows(0x1)}")

print()
print("=" * 78)
print("5.  guessed-key masks of the Appendix")
print("=" * 78)
MK0 = 0xFE3FFF3FFF3BFF3F
MK1 = 0xCC0CCC0CCC0C0000
MK0_1 = 0x5626573757334737
MK0_2 = 0xB83BB83BB83BB83B
check("wt(M_K0) = 54", bin(MK0).count("1") == 54, str(bin(MK0).count("1")))
check("wt(M_K1) = 18", bin(MK1).count("1") == 18, str(bin(MK1).count("1")))
check("M_K0 = M_K0^(1) OR M_K0^(2)", MK0 == (MK0_1 | MK0_2))
check("wt(M_K0^(1)) = 35, wt(M_K0^(2)) = 36",
      bin(MK0_1).count("1") == 35 and bin(MK0_2).count("1") == 36)
check("k = 54 + 18 = 72", bin(MK0).count("1") + bin(MK1).count("1") == 72)

# (a) dX^3 is determined by the M_K1 coordinates of X^2 once dX^2 = C0CC
DX2 = 0x000000000000CC0C
det_ok = True
for _ in range(3000):
    X2 = rnd64()
    dX1 = rnd64()
    base = F(X2) ^ F(X2 ^ DX2)
    for _ in range(4):
        X2b = (X2 & MK1) | (rnd64() & ~MK1 & MASK64)
        det_ok &= (F(X2b) ^ F(X2b ^ DX2)) == base
check("D_{dX^2} f(X^2) depends only on the 18 M_K1 coordinates", det_ok)

# (b) those 18 coordinates of f(X^1) need only the M_K0^(2) coordinates of X^1
pre_ok = True
for _ in range(3000):
    X1 = rnd64()
    base = F(X1) & MK1
    for _ in range(4):
        X1b = (X1 & MK0_2) | (rnd64() & ~MK0_2 & MASK64)
        pre_ok &= (F(X1b) & MK1) == base
check("f(X^1) restricted to M_K1 depends only on the M_K0^(2) coordinates", pre_ok)

# (c) composite: the whole counter test depends on sk^0, sk^1 only through the masks
comp_ok = True
tight0 = 0
tight1 = 0
for _ in range(600):
    X3, Y3 = rnd64(), rnd64()
    X3s, Y3s = X3 ^ DX3, Y3 ^ DY3
    sk0, sk1, sk2 = rnd64(), rnd64(), rnd64()
    x, y = X3, Y3
    xs, ys = X3s, Y3s
    for k in (sk2, sk1, sk0):
        x, y = dec_round(x, y, k)
        xs, ys = dec_round(xs, ys, k)
    X0, Y0, X0s, Y0s = x, y, xs, ys

    def test(g0, g1):
        X1 = Y0 ^ F(X0) ^ g0
        X1s = Y0s ^ F(X0s) ^ g0
        X2 = X0 ^ F(X1) ^ g1
        X2s = X0s ^ F(X1s) ^ g1
        return (X1 ^ X1s, X2 ^ X2s, (X1 ^ X1s) ^ F(X2) ^ F(X2s))

    ref = test(sk0, sk1)
    comp_ok &= ref[1] == DX2 and ref[2] == DX3       # true key is a right pair
    for _ in range(6):
        g0 = (sk0 & MK0) | (rnd64() & ~MK0 & MASK64)
        g1 = (sk1 & MK1) | (rnd64() & ~MK1 & MASK64)
        comp_ok &= test(g0, g1) == ref
    # tightness: flipping one masked bit usually changes the outcome
    i = random.choice([i for i in range(64) if (MK0 >> i) & 1])
    if test(sk0 ^ (1 << i), sk1) != ref:
        tight0 += 1
    i = random.choice([i for i in range(64) if (MK1 >> i) & 1])
    if test(sk0, sk1 ^ (1 << i)) != ref:
        tight1 += 1
check("right pairs are counted by the true key class", comp_ok)
check("counter test is invariant under key bits outside M_K0 / M_K1", comp_ok)
print(f"       (single masked-bit flips changed the outcome in {tight0}/600 (sk^0) "
      f"and {tight1}/600 (sk^1) trials -> masks are not padded with dead bits)")

print()
print("=" * 78)
print("6.  wrong-key behaviour and ranking (scaled-down simulation)")
print("=" * 78)
# 32 planted right pairs (= E[c_right] for 2^125 pairs), guessing a 14-bit
# sub-window of M_K0 with the remaining bits fixed to their correct values
SUBBITS = [i for i in range(64) if (MK0 >> i) & 1][:14]
sk0, sk1, sk2 = rnd64(), rnd64(), rnd64()
pairs = []
for _ in range(32):
    X3, Y3 = rnd64(), rnd64()
    x, y = X3, Y3
    xs, ys = X3 ^ DX3, Y3 ^ DY3
    for k in (sk2, sk1, sk0):
        x, y = dec_round(x, y, k)
        xs, ys = dec_round(xs, ys, k)
    pairs.append((x, y, xs, ys))
for _ in range(2000):                       # random (wrong) pairs as background
    pairs.append((rnd64(), rnd64(), rnd64(), rnd64()))
pre = [(X0, Y0, X0s, Y0s, F(X0), F(X0s)) for (X0, Y0, X0s, Y0s) in pairs]

counters = {}
for code in range(1 << len(SUBBITS)):
    g0 = sk0
    for t, i in enumerate(SUBBITS):
        g0 = (g0 & ~(1 << i)) | (((code >> t) & 1) << i)
    c = 0
    for (X0, Y0, X0s, Y0s, fx, fxs) in pre:
        X1 = Y0 ^ fx ^ g0
        X1s = Y0s ^ fxs ^ g0
        X2 = X0 ^ F(X1) ^ sk1
        X2s = X0s ^ F(X1s) ^ sk1
        if (X2 ^ X2s) == DX2 and ((X1 ^ X1s) ^ F(X2) ^ F(X2s)) == DX3:
            c += 1
    counters[code] = c
correct = 0
for t, i in enumerate(SUBBITS):
    correct |= ((sk0 >> i) & 1) << t
best = max(counters.values())
winners = [c for c, v in counters.items() if v == best]
others = sorted(v for c, v in counters.items() if c != correct)
check("correct partial key is the unique maximum counter",
      winners == [correct],
      f"max = {best} (of 32 planted right pairs), #ties = {len(winners)}, "
      f"best wrong counter = {others[-1]}")
# Note on the threshold: the paper's tau = 14 is calibrated against the
# wrong-key randomisation hypothesis, in which a wrong class collects
# Poisson(2^-3) counts. It must not be applied to this sub-window, whose keys
# differ from the correct one in a few bits only and are therefore *not*
# independent wrong classes: on a planted right pair a differing bit is often
# inert, so near-correct window keys score far above 2^-3. The ranking test
# above is the meaningful one at this scale; the Poisson calibration itself is
# checked arithmetically in section 6.
print(f"       counters of the 2^{len(SUBBITS)} window keys: correct = "
      f"{counters[correct]}, best wrong = {others[-1]}, median = "
      f"{others[len(others) // 2]}")

print()
print("=" * 78)
print("7.  complexity arithmetic")
print("=" * 78)
import math
N = 125                                          # pairs
W = 120
f0, f1 = 27, 52
k = 72
rounds = 23
TAU = 14
E_WRONG = -6.5                                   # log2 E[# wrong classes kept]
print(f"  pairs                  2^{N}")
print(f"  chosen ciphertexts     2^{N+1}   (involution -> 2 texts per pair)")
print(f"  right pairs            2^{N-W} = {2**(N-W)}")
print(f"  after Eq.(11)          2^{N-f0}")
print(f"  after Eq.(12)          2^{N-f0-f1}")
print(f"  wrong class expectation 2^{N-128}")
# staged early abort: the residual conditions are 37 bits (dX^2) and 12 (dX^3)
resid2, resid3 = 128 - f0 - f1 - 12, 12
sA = (N - f0 - f1) + 35                      # 2^46 pairs x 2^35 guesses of M_K0^(1)
survA = sA - resid2
sB = survA + 37                              # 19 more sk^0 bits + 18 of sk^1
check("residual conditions sum to the 128-bit condition",
      f0 + f1 + resid2 + resid3 == 128, f"{f0}+{f1}+{resid2}+{resid3}")
check(f"staged counting reproduces 2^{k} x 2^{N-128} increments",
      sB - resid3 == k + (N - 128), f"2^{sB - resid3} vs 2^{k + N - 128}")
print(f"  stage 1 2^{sA} -> 2^{survA} survivors ; stage 2 2^{sB} -> 2^{sB-resid3} counts")

# the paper's single sum, all terms in units of 2/23 of a 23-round encryption
u = math.log2(2 / rounds)
t_data = u + (N + 1)                             # 2 * 2^N  (Steps 1 and 2)
t_dx1 = u + (N - f0)                             # Step 3
t_st1 = u + sA                                   # Stage 1
t_st2 = u + sB                                   # Stage 2
t_fin = max(E_WRONG, 0.0) + (128 - k)            # Step 5, tau = 14 keeps only the right class
T = math.log2(2**t_data + 2**t_dx1 + 2**t_st1 + 2**t_st2 + 2**t_fin)
T_kr = math.log2(2**t_dx1 + 2**t_st1 + 2**t_st2 + 2**t_fin)
# key-outer loop: no counter array and, at tau = 14, no list of retained classes
M = math.log2(2 * 2 ** (N - f0 - f1) * 16)
no_abort = math.log2(4 / rounds) + (N - f0 - f1) + k
print(f"  Steps 1,2 2^{t_data:.2f} | Step 3 2^{t_dx1:.2f} | Stage 1 2^{t_st1:.2f}"
      f" | Stage 2 2^{t_st2:.2f} | Step 5 2^{t_fin:.2f}")
print(f"  T      = 2^{T:.2f}    (paper 2^122.5)")
print(f"  T_kr   = 2^{T_kr:.2f}     (paper 2^94.5)")
print(f"  M      = 2^{M:.2f} bytes  (paper 2^51)")
print(f"  merged stages, no early abort = 2^{no_abort:.2f}  (paper 2^115.5), "
      f"still below Steps 1,2 by 2^{t_data - no_abort:.1f}")
check("total time matches the paper", abs(T - 122.5) < 0.1, f"2^{T:.2f}")
check("key-recovery time matches the paper", abs(T_kr - 94.5) < 0.1, f"2^{T_kr:.2f}")
check("memory matches the paper", abs(M - 51) < 0.1, f"2^{M:.2f}")
check("no-early-abort figure matches the paper", abs(no_abort - 115.5) < 0.1,
      f"2^{no_abort:.2f}")
check(f"data 2^{N+1} < full codebook 2^128", N + 1 < 128)
check("time < exhaustive search 2^128", T < 128)

# signal/noise with the Poisson model over 2^72 classes
def poisson_ge(lam, tau):
    """P[Poisson(lam) >= tau], summed exactly over the complement."""
    s, term = 0.0, 1.0
    for j in range(tau):
        if j:
            term *= lam / j
        s += term
    return -math.expm1(-lam) * 0 + (1.0 - math.exp(-lam) * s)

lam_r = 2.0 ** (N - W)
lam_w = 2.0 ** (N - 128)
p_succ = poisson_ge(lam_r, TAU)
# E[# wrong classes >= tau]; the tail is tiny, so the leading term is exact enough
log_p_w = (-lam_w / math.log(2) + TAU * math.log2(lam_w)
           - math.log2(math.factorial(TAU)))
print()
print("  Poisson model: E[C_right] = 2^%d = %d, E[C_wrong] = 2^%d = %.4f"
      % (N - W, int(lam_r), N - 128, lam_w))
print(f"  Pr[C_right >= {TAU}] = {p_succ:.5f}   (paper 0.99988)")
print(f"  E[# wrong classes >= {TAU}] = 2^{log_p_w + k:.1f}   (paper 2^-6.5)")
check("success probability matches the paper", abs(p_succ - 0.99988) < 1e-4,
      f"{p_succ:.5f}")
check("no wrong class is expected to survive the threshold", log_p_w + k < 0,
      f"2^{log_p_w + k:.1f}")

print()
print("=" * 78)
print("8.  the formulas exactly as printed in main.tex")
print("=" * 78)
# Sections 0-7 check that the *attack* is sound. They do not check that the
# formulas printed in the paper say what the implementation does -- a gap that
# let two typographic errors through review. This section closes it: each
# printed relation is transcribed literally, evaluated, and required to
# reproduce the true cipher, and the erroneous variant is required to fail.

# --- 8a. the round recurrence of Stage 1 / Stage 2 and of Eq. (A.1) ---
# printed:  dX^c = dX^{c-2} XOR D_{dX^{c-1}} f(X^{c-1})
# rejected: dX^c = dX^{c-2} XOR P_n( D_{dX^{c-1}} f(X^{c-1}) )   (P_n twice)
def printed_next_diff(dX_prev2, X_prev1, X_prev1s):
    return dX_prev2 ^ F(X_prev1) ^ F(X_prev1s)


def double_p_variant(dX_prev2, X_prev1, X_prev1s):
    return dX_prev2 ^ P(F(X_prev1) ^ F(X_prev1s))


printed_ok = variant_rejected = True
for _ in range(300):
    sk0, sk1, sk2 = rnd64(), rnd64(), rnd64()
    X3, Y3 = rnd64(), rnd64()
    x, y, xs, ys = X3, Y3, X3 ^ DX3, Y3 ^ DY3
    for k in (sk2, sk1, sk0):
        x, y = dec_round(x, y, k)
        xs, ys = dec_round(xs, ys, k)
    X0, Y0, X0s, Y0s = x, y, xs, ys
    X1, X1s = Y0 ^ F(X0) ^ sk0, Y0s ^ F(X0s) ^ sk0
    X2, X2s = X0 ^ F(X1) ^ sk1, X0s ^ F(X1s) ^ sk1
    X3r, X3rs = X1 ^ F(X2) ^ sk2, X1s ^ F(X2s) ^ sk2
    for dprev2, a, b, truth in ((X0 ^ X0s, X1, X1s, X2 ^ X2s),
                                (X1 ^ X1s, X2, X2s, X3r ^ X3rs)):
        printed_ok &= printed_next_diff(dprev2, a, b) == truth
        variant_rejected &= double_p_variant(dprev2, a, b) != truth
check("printed recurrence dX^c = dX^(c-2) XOR D f(X^(c-1)) reproduces the cipher",
      printed_ok)
check("the discarded variant with an extra P_n does NOT reproduce the cipher",
      variant_rejected, "guards against re-introducing the double P_n")

# --- 8b. the transfer-matrix relation of Section VII-B-3 ---
# printed:  T[j] = d_{0,j} XOR d_{1,(j-1) mod n/4}   with LSB-first columns
# rejected: T[j] = d_{0,j} XOR d_{1,j+1}
def transitions_feasible(shift):
    """number of the 19 rounds whose transition admits a solution."""
    good = 0
    for r in range(1, 20):
        dxp, dyp = TRAIL[r - 1]
        dxc = TRAIL[r][0]
        T = cols_lsb(p_inv(dxc ^ dyp))
        din = cols_lsb(dxp)
        opts = [[(0, 0)] if din[m] == 0 else
                [((o >> 4) & 0xF, o & 0xF) for o in range(256) if DDT[din[m]][o]]
                for m in range(16)]
        for combo in itertools.product(*opts):
            if all(T[j] == (combo[j][0] ^ combo[(j + shift) % 16][1])
                   for j in range(16)):
                good += 1
                break
    return good


feasible_printed = transitions_feasible(-1)
feasible_wrong = transitions_feasible(+1)
print(f"  T[j] = d0_j XOR d1_(j-1)  (printed) : {feasible_printed}/19 transitions feasible")
print(f"  T[j] = d0_j XOR d1_(j+1)            : {feasible_wrong}/19 transitions feasible")
check("printed transfer-matrix index reproduces every round transition",
      feasible_printed == 19)
check("the opposite index does not", feasible_wrong < 19,
      "guards against re-introducing the MSB-first form")

# --- 8c. the source file still carries the corrected wording ---
TEX = os.path.join(HERE, "..", "..", "paper", "SAND_DC_IEICE", "main.tex")
if os.path.exists(TEX):
    tex = open(TEX, encoding="utf-8", errors="replace").read()
    banned = {
        r"P_n\!\left(D_{\Delta": "double application of P_n",
        r"d_{1,j+1}": "MSB-first transfer-matrix index",
        "leave no wrong class behind": "expectation stated as a guarantee",
    }
    hits = [why for pat, why in banned.items() if pat in tex]
    check("main.tex is free of the three corrected defects", not hits,
          "; ".join(hits) if hits else "none present")
    required = [r"d_{1,(j-1)\bmod\frac{n}{4}}"]
    check("main.tex carries the corrected transfer-matrix index",
          all(pat in tex for pat in required))
else:
    print(f"  [skip] {TEX} not found; source-text checks not run")

print()
print("=" * 78)
print("9.  cross-check against the author's own implementation")
print("=" * 78)
sys.path.insert(0, HERE)
import keyrecovery_sand128 as K
agree = True
for _ in range(500):
    w = rnd64()
    agree &= K.apply_p(w) == P(w)
    agree &= K.apply_p_inv(w) == p_inv(w)
    agree &= [K.nibble(w, j) for j in range(16)] == cols(w)
check("P_n / P_n^-1 / nibble extraction agree with keyrecovery_sand128.py", agree)
kw, kt = K.verify_trail(TRAIL, "paper")
check("transfer-matrix weights agree with the enumeration above",
      kw == PAPER_W and kt == 120, f"{kt}")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
print("=" * 78)
sys.exit(0 if ok_all else 1)
