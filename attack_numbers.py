#!/usr/bin/env python3
"""
論文に載せる攻撃の数値をすべて確定させる。

構成 (論文 SAND-64 の 16 段攻撃を SAND-128 に移したもの):
    第 0..2 段     : 平文側に付ける 3 段 (鍵推測)
    第 3..22 段    : 19 段差分識別子 (重み 120)
    第 22..23 段   : 暗号文側 1 段 (最終ラウンド鍵に依存しない対応付け)
    合計 23 段
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keyrecovery_sand128 as K
import filters_sand128 as F
import keybits_sand128 as KB

R = 19
R_IN = 3
R_OUT = 1
W = 120
MARGIN = 4
BLOCK = 128
KEYLEN = 128
NCOLS = 16


def load_trails(path):
    trails, cur = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith("#"):
            if cur:
                trails.append(cur)
            cur = []
        elif line:
            p = line.split()
            cur.append((int(p[1], 16), int(p[2], 16)))
    if cur:
        trails.append(cur)
    return [t for t in trails if len(t) == R + 1]


def evaluate(trail):
    dx_s = trail[R][1]          # ΔY^R -> 攻撃の ΔX^{R_IN}
    dy_s = trail[R][0]          # ΔX^R -> 攻撃の ΔY^{R_IN}
    seq = F.truncated_backward(dx_s, dy_s, R_IN, NCOLS, K.P_BOX)
    dx0, _ = seq[R_IN]
    dx1, _ = seq[R_IN - 1]
    m0, v0 = F.known_mask(dx0, NCOLS)
    m1, v1 = F.known_mask(dx1, NCOLS)
    f_bits = bin(m0).count("1") + bin(m1).count("1")
    k, masks = KB.key_bits(dx_s, dy_s, R_IN, NCOLS, K.P_BOX)
    return dict(dx_s=dx_s, dy_s=dy_s, m0=m0, v0=v0, m1=m1, v1=v1,
                F=f_bits, k=k, masks=masks)


def main():
    trails = load_trails(os.path.join(HERE, "trails_sand128_r19_w120.txt"))
    print(f"{len(trails)} verified 19-round characteristics of weight {W}")

    scored = []
    for t in trails:
        e = evaluate(t)
        # 生き残り対 2^(W+MARGIN-F) を最小に, 次いで k を小さく
        scored.append((-e["F"], e["k"], t, e))
    scored.sort(key=lambda z: (z[0], z[1]))
    _, _, trail, e = scored[0]

    exact_w, total = K.verify_trail(trail, "chosen")
    assert total == W

    print("\n===== selected characteristic =====")
    print(f"{'r':>3} {'dX^r':>17} {'dY^r':>17} {'log2 DP':>8} {'log2 DCP':>9}")
    acc = 0
    for r, (dx, dy) in enumerate(trail):
        if r == 0:
            print(f"{r:>3} {K.fmt(dx):>17} {K.fmt(dy):>17} {'':>8} {'':>9}")
        else:
            acc += exact_w[r - 1]
            print(f"{r:>3} {K.fmt(dx):>17} {K.fmt(dy):>17} {-exact_w[r-1]:>8} {-acc:>9}")
    print(f"exact weight (transfer matrix) = {total}  ->  p = 2^-{total}")

    ins = K.insert_backward(trail, R_IN)
    ct = ins[R_IN + R]
    pt = ins[R_IN]
    total_rounds = R_IN + R + R_OUT
    print(f"\n===== insertion (chosen ciphertext), {total_rounds}-round attack =====")
    print(f"(dX^{R_IN+R}, dY^{R_IN+R}) -> (dX^{R_IN}, dY^{R_IN})")
    print(f"(dX^{R_IN+R}, dY^{R_IN+R}) = ({K.fmt(ct[0])}, {K.fmt(ct[1])})   <- a, b")
    print(f"(dX^{R_IN}, dY^{R_IN}) = ({K.fmt(pt[0])}, {K.fmt(pt[1])})")

    print(f"\n===== key-free filters =====")
    print(f"dX^0 & {K.fmt(e['m0'])} = {K.fmt(e['v0'])}   wt = {bin(e['m0']).count('1')}")
    print(f"dX^1 & {K.fmt(e['m1'])} = {K.fmt(e['v1'])}   wt = {bin(e['m1']).count('1')}")
    print(f"F = {e['F']}")

    print(f"\n===== round-key masks =====")
    for t in sorted(e["masks"]):
        m = e["masks"][t]
        print(f"M_K{t-1} (sk^{t-1}) = {K.fmt(m)}   wt = {bin(m).count('1')}")
    print(f"k = {e['k']}")

    # ---- complexities ----
    N = W + MARGIN                      # log2 対数
    D = N + 1                           # 選択暗号文数
    f0 = bin(e["m0"]).count("1")
    f1 = bin(e["m1"]).count("1")
    surv1 = N - f0
    surv2 = surv1 - f1
    k = e["k"]
    rounds = total_rounds

    import math
    t1 = math.log2(2 * (1 / rounds)) + surv1
    t2 = math.log2(2 * (2 / rounds)) + surv2 + k
    t3 = KEYLEN - k
    T = max(t1, t2, t3) + math.log2(
        2 ** (t1 - max(t1, t2, t3)) + 2 ** (t2 - max(t1, t2, t3)) + 2 ** (t3 - max(t1, t2, t3)))
    mem = math.log2(2 ** k + 2 * 2 ** surv2 * BLOCK / 8)

    print(f"\n===== complexities ({rounds}-round attack) =====")
    print(f"pairs N            = 2^{N}")
    print(f"data D             = 2^{D} chosen ciphertexts")
    print(f"right pairs        = 2^{N - W}")
    print(f"after filter 1     = 2^{surv1}")
    print(f"after filter 2     = 2^{surv2}")
    print(f"noise per class    = 2^{N - BLOCK}  (full {BLOCK}-bit condition)")
    print(f"time term 1        = 2^{t1:.2f}")
    print(f"time term 2        = 2^{t2:.2f}")
    print(f"time term 3        = 2^{t3}   (exhaustive over {KEYLEN}-{k} bits)")
    print(f"time T             = 2^{T:.2f}")
    print(f"memory M           = 2^{mem:.2f} bytes")


if __name__ == "__main__":
    main()
