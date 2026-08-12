#!/usr/bin/env python3
"""Check the quantitative claims of SAND_review_memo.txt."""
import os, sys, math, itertools, random, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keyrecovery_sand128 as K
import filters_sand128 as FL
import keybits_sand128 as KB
import attack_numbers as A
# --- primitives inlined (same code as verify_sand128_key_recovery.py) ---
P16 = [14, 15, 8, 9, 2, 3, 12, 13, 6, 7, 0, 1, 10, 11, 4, 5]
MASK64 = (1 << 64) - 1
def G0(u):
    u3,u2,u1,u0 = (u>>3)&1,(u>>2)&1,(u>>1)&1,u&1
    y0=(u3&u2)^u0; y3=(y0&u1)^u3
    return (y3<<3)|(u2<<2)|(u1<<1)|y0
def G1(u):
    u3,u2,u1,u0 = (u>>3)&1,(u>>2)&1,(u>>1)&1,u&1
    y2=(u3&u1)^u2; y1=(y2&u0)^u1
    return (u3<<3)|(y2<<2)|(y1<<1)|u0
DDT = [[0]*256 for _ in range(16)]
for _d in range(16):
    for _u in range(16):
        DDT[_d][((G0(_u)^G0(_u^_d))<<4)|(G1(_u)^G1(_u^_d))] += 1
PTAB = [0]*65536
for _v in range(65536):
    _o = 0
    for _j in range(16):
        if (_v >> (15-_j)) & 1:
            _o |= 1 << P16[15-_j]
    PTAB[_v] = _o
def F(x):
    r3,r2,r1,r0 = (x>>48)&0xFFFF,(x>>32)&0xFFFF,(x>>16)&0xFFFF,x&0xFFFF
    a0=(r3&r2)^r0; a1=r1; a2=r2; a3=(a0&r1)^r3
    s3=((r3<<1)|(r3>>15))&0xFFFF; s2=((r2<<1)|(r2>>15))&0xFFFF
    s1=((r1<<1)|(r1>>15))&0xFFFF; s0=((r0<<1)|(r0>>15))&0xFFFF
    b3=s3; b2=(s3&s1)^s2; b1=(b2&s0)^s1; b0=s0
    return (PTAB[a3^b3]<<48)|(PTAB[a2^b2]<<32)|(PTAB[a1^b1]<<16)|PTAB[a0^b0]
def dec_round(xn, yn, sk):
    x = yn
    return x & MASK64, (xn ^ F(x) ^ sk) & MASK64
def hx(w): return f"{w:016X}"

random.seed(7)
rnd = lambda: random.getrandbits(64)
DX3, DY3 = 0x0000000000000510, 0x000000000000C0CC
DX2 = 0x000000000000C0CC
M0, V0 = 0xE222C546C444F333, 0x000000000000F333
M1, V1 = 0xCCCFFFFF9D9BFFFF, 0x0000000000005444

print("=" * 78)
print("[1] SSb weight distribution and the marginal-vs-joint example")
print("=" * 78)
wts = collections.Counter()
for a in range(1, 16):
    for o in range(256):
        if DDT[a][o]:
            wts[4 - int(math.log2(DDT[a][o]))] += 1
print("  nonzero-input transition weights:", dict(sorted(wts.items())),
      "-> memo says {2,3} for active S-boxes (plus weight 0 for the inactive one)")
mg0 = collections.Counter(); mg1 = collections.Counter()
for u in range(16):
    mg0[G0(u) ^ G0(u ^ 0x6)] += 1
    mg1[G1(u) ^ G1(u ^ 0x6)] += 1
print(f"  marginal DDT_N0[0x6][0x6] = {mg0[0x6]}, DDT_N1[0x6][0x2] = {mg1[0x2]}, "
      f"joint DDT_SSb[0x6][0x62] = {DDT[0x6][0x62]}  (memo: 4, 8, 0)")

print()
print("=" * 78)
print("[2] B[r] versus the designers' active-S-box counts (SAND Table 10)")
print("=" * 78)
designers = {2:1,3:2,4:4,5:6,6:9,7:12,8:16,9:20,10:25,11:30,12:34}
import json
B = {int(k): v for k, v in
     json.load(open(os.path.join(HERE, "bounds_sand128.json"))).items()}
ok = all(B[r] == 2 * designers[r] for r in designers)
print("   r      :", "".join(f"{r:>5}" for r in sorted(designers)))
print("   2 x AS :", "".join(f"{2*designers[r]:>5}" for r in sorted(designers)))
print("   B[r]   :", "".join(f"{B[r]:>5}" for r in sorted(designers)))
print(f"  -> {'MATCH (memo confirmed)' if ok else 'MISMATCH'}")

print()
print("=" * 78)
print("[3] memo [F]: superadditivity  B[22] >= B[20] + B[2] >= 128")
print("=" * 78)
print(f"  B[20] >= 126 (w<=125 UNSAT in sand128_bounds.log), B[2] = {B[2]}")
print(f"  => B[22] >= 128, i.e. the 2^-128 bound is reached at 22 rounds,")
print(f"     one round earlier than the designers' estimate of 23 rounds")
print(f"     (SAND paper: 11-round 30 AS + 12-round 34 AS = 64 AS = weight 128).")
print(f"  best bound for 21 rounds: B[21] >= B[20]+B[1] = {126+B[1]}  (< 128, so 22 is tight here)")

print()
print("=" * 78)
print("[4] memo [A]: residual condition after the two key-free filters (36 + 12?)")
print("=" * 78)
# build a right pair so that dX^0, dX^1 obey the filters
def right_pair():
    while True:
        X3, Y3 = rnd(), rnd()
        sk = [rnd(), rnd(), rnd()]
        x, y, xs, ys = X3, Y3, X3 ^ DX3, Y3 ^ DY3
        for r in (2, 1, 0):
            x, y = dec_round(x, y, sk[r]); xs, ys = dec_round(xs, ys, sk[r])
        return x, y, xs, ys, sk

X0, Y0, X0s, Y0s, sk = right_pair()
dX0, dY0 = X0 ^ X0s, Y0 ^ Y0s
dX1 = dY0 ^ F(X0) ^ F(X0s)
assert (dX0 & M0) == V0 and (dX1 & M1) == V1

# how many bits of dX^2 are already fixed once the pair passed the filters?
and2, or2 = MASK64, 0
for _ in range(4000):
    g0 = rnd()
    X1, X1s = Y0 ^ F(X0) ^ g0, Y0s ^ F(X0s) ^ g0
    d2 = dX0 ^ F(X1) ^ F(X1s)
    and2 &= d2; or2 |= d2
const2 = ~(and2 ^ or2) & MASK64
print(f"  bits of dX^2 constant over all sk^0 guesses : {bin(const2).count('1')}"
      f"   (memo: 28 already consumed by Eq.(11))")
print(f"  -> residual condition on dX^2 = 64 - {bin(const2).count('1')} = "
      f"{64 - bin(const2).count('1')} bits   (memo says 36)")
print(f"  those constant bits agree with C0CC : {(and2 & const2) == (DX2 & const2)}")

and3, or3 = MASK64, 0
for _ in range(4000):
    X2 = rnd()
    d3 = dX1 ^ F(X2) ^ F(X2 ^ DX2)
    and3 &= d3; or3 |= d3
const3 = ~(and3 ^ or3) & MASK64
print(f"  bits of dX^3 constant over all X^2 (given dX^2=C0CC) : {bin(const3).count('1')}"
      f"   (memo: 52 consumed by Eq.(12))")
print(f"  -> residual condition on dX^3 = 64 - {bin(const3).count('1')} = "
      f"{64 - bin(const3).count('1')} bits   (memo says 12)")
print(f"  those constant bits agree with 0510 : {(and3 & const3) == (DX3 & const3)}")
tot = (64 - bin(const2).count('1')) + (64 - bin(const3).count('1'))
print(f"  total residual = {tot} bits   (must be 128 - 28 - 52 = 48)")

print()
print("  staged early-abort cost:")
sA = 44 + 37
print(f"    stage A : 2^44 pairs x 2^37 (M_K0^(1))            = 2^{sA}"
      f"  -> 2^{sA} x 2^-{64-bin(const2).count('1')} = 2^{sA-(64-bin(const2).count('1'))} survivors")
sB = (sA - (64 - bin(const2).count('1'))) + 40
print(f"    stage B : 2^{sA-(64-bin(const2).count('1'))} x 2^40 (22 more sk^0 bits + M_K1) = 2^{sB}"
      f"  -> counts 2^{sB-(64-bin(const3).count('1'))}")
print(f"    total counts 2^{sB-(64-bin(const3).count('1'))} vs 2^77 classes x 2^-4 = 2^73"
      f"  -> {'consistent' if sB-(64-bin(const3).count('1'))==73 else 'INCONSISTENT'}")
print(f"    in encryptions: 2^{sB} x 4/23 = 2^{sB + math.log2(4/23):.1f}  (memo: 2^82.5)")

print()
print("=" * 78)
print("[5] cost that CANNOT be avoided: scanning the 2^124 pairs")
print("=" * 78)
ops_per_enc = 23 * 15          # ~15 64-bit word operations per round
print(f"  one 23-round encryption ~ {ops_per_enc} word ops = 2^{math.log2(ops_per_enc):.1f}")
for cost, lbl in ((3, "XOR+AND+compare"), (15, "one full round"), (ops_per_enc, "one encryption")):
    print(f"    filtering 2^124 pairs at {lbl:>16} each : "
          f"2^{124 + math.log2(cost/ops_per_enc):.1f} encryptions")
print("  => after the early abort of [A], THIS is the dominant term, not 2^92.5.")
print(f"  memo [B] claims offline = 2^92.5; that ignores the 2^124 scan.")

print()
print("=" * 78)
print("[6] memo [D]: differential clustering among the 151 weight-120 trails")
print("=" * 78)
trails = A.load_trails(os.path.join(HERE, "trails_sand128_r19_w120.txt"))
groups = collections.Counter((t[0][0], t[0][1], t[19][0], t[19][1]) for t in trails)
print(f"  {len(trails)} trails -> {len(groups)} distinct (input,output) differentials")
print(f"  largest cluster: {max(groups.values())} trails sharing one differential")
sel = A.load_trails(os.path.join(HERE, "trails_sand128_r19_w120.txt"))
e = [t for t in trails if (t[19][0], t[19][1]) == (0x000000000000C0CC, 0x0000000000000510)
     and (t[0][0], t[0][1]) == (0x0000000000001115, 0x000000000000CC0C)]
print(f"  trails sharing the endpoints of the selected characteristic: {len(e)}")
for k, v in sorted(groups.items(), key=lambda z: -z[1])[:5]:
    print(f"    {hx(k[0])} {hx(k[1])} -> {hx(k[2])} {hx(k[3])} : {v}")

print()
print("=" * 78)
print("[7] memo [C]: can a 4th round be prepended?")
print("=" * 78)
for r_in in (3, 4):
    seq = FL.truncated_backward(DX3, DY3, r_in, 16, K.P_BOX)
    known = [bin(FL.known_mask(seq[t][0], 16)[0]).count("1") for t in range(r_in + 1)]
    k, masks = KB.key_bits(DX3, DY3, r_in, 16, K.P_BOX)
    print(f"  prepending {r_in} rounds: known bits of dX^(s-t) = {known[::-1]}"
          f"  -> total filter {sum(known[1:])} bits, guessed key k = {k}")
print("  (k is capped by the 128-bit master key: sk^0||sk^1 = K^1||K^0.)")
