#!/usr/bin/env python3
"""
論文誌版 (paper/SAND_DC_IEICE) の SAND-128 23 段鍵回復攻撃の数値を,
特性データから導出して確定させる。

構成:
    第 0..2 段   : 平文側に付ける 3 段。sk^0 の 54 ビットと sk^1 の 18 ビットを推測する。
    第 3..22 段  : 19 段差分識別子 (重み 120) を **前向きに** 挿入する。すなわち
                   攻撃の (ΔX^3, ΔY^3) が特性の入力 (ΔX^0, ΔY^0) に,
                   (ΔX^22, ΔY^22) が特性の出力 (a, b) に一致する。
    第 22..23 段 : 暗号文側 1 段。(ΔX^22, ΔY^22) = (a, b) を強制する暗号文の
                   対応付けは最終ラウンド鍵に依存しない不動点自由対合なので無料。
    合計 23 段

注意: 識別子は前向きに挿入する。特性を反転して挿入する
keyrecovery_sand128.insert_backward() を使うと, 得られる鍵非依存フィルタは
右対の 100% が違反する成立しない条件になる。

導出した数値は verify_sand128_key_recovery.py (仕様から書き起こした独立実装)
および論文の Table 13 と一致することを最後に検査する。
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
import keyrecovery_sand128 as K
import filters_sand128 as F
import keybits_sand128 as KB

R = 19                  # 識別子の段数
R_IN = 3                # 平文側に付ける段数
R_OUT = 1               # 暗号文側に付ける段数
W = 120                 # 識別子の重み -log2 p
N_PAIRS = 125           # log2 の対数。E[C_right] = 2^5 = 32 となり閾値 tau=14 が効く
BLOCK = 128
KEYLEN = 128
NCOLS = 16
TAU = 14                # カウンタの閾値 (論文 Section V)
STAGE1_BITS = 35        # wt(M_K0^(1)): 第 1 段階で推測する sk^0 のビット数 (Appendix)
RESID3 = 12             # ΔX^3 に残る条件のビット数 (第 2 段階で消費する)
LOG_E_WRONG = -6.5      # log2 E[閾値を越える誤鍵クラス数]

# 論文 Table 13 および Appendix の公表値
PAPER = dict(f0=27, f1=52, k=72, M_K0=54, M_K1=18,
             data=126, T=122.5, T_kr=94.5, mem=51.0, p_succ=0.99988)
PAPER_TRAIL_IN = (0x0000000000001115, 0x000000000000CC0C)
PAPER_TRAIL_OUT = (0x000000000000C0CC, 0x0000000000000510)


def wt(x):
    return bin(x).count("1")


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
    """前向き挿入したときの鍵非依存フィルタと推測鍵ビット数を求める。"""
    dx_s, dy_s = trail[0]           # 攻撃の (ΔX^{R_IN}, ΔY^{R_IN})
    seq = F.truncated_backward(dx_s, dy_s, R_IN, NCOLS, K.P_BOX)
    m0, v0 = F.known_mask(seq[R_IN][0], NCOLS)          # ΔX^0 上の条件
    m1, v1 = F.known_mask(seq[R_IN - 1][0], NCOLS)      # ΔX^1 上の条件
    k, masks = KB.key_bits(dx_s, dy_s, R_IN, NCOLS, K.P_BOX)
    return dict(dx_s=dx_s, dy_s=dy_s, m0=m0, v0=v0, m1=m1, v1=v1,
                F=wt(m0) + wt(m1), k=k, masks=masks)


def poisson_ge(lam, tau):
    """P[Poisson(lam) >= tau]"""
    s, term = 0.0, 1.0
    for j in range(1, tau):
        term *= lam / j
        s += term
    return 1.0 - math.exp(-lam) * (1.0 + s)


ok_all = True


def check(name, cond, extra=""):
    global ok_all
    ok_all &= bool(cond)
    print(f"[{'OK ' if cond else 'FAIL'}] {name}{(' : ' + extra) if extra else ''}")
    return cond


def main():
    trails = load_trails(os.path.join(RESULTS, "trails_sand128_r19_w120.txt"))
    print(f"{len(trails)} verified 19-round characteristics of weight {W}")

    # 生き残り対 2^(W+margin-F) を最小に, 次いで推測鍵ビット数 k を小さく
    scored = sorted((-e["F"], e["k"], i, e)
                    for i, e in ((i, evaluate(t)) for i, t in enumerate(trails)))
    _, _, idx, e = scored[0]
    trail = trails[idx]

    best = {(-a, b) for a, b, _, _ in scored[:1]}
    n_best = sum(1 for a, b, _, _ in scored if (-a, b) in best)
    print(f"best (F, k) = ({e['F']}, {e['k']}), attained by {n_best} of them")

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

    total_rounds = R_IN + R + R_OUT
    a, b = trail[R]
    print(f"\n===== insertion (forward), {total_rounds}-round attack =====")
    print(f"(dX^{R_IN}, dY^{R_IN})   = ({K.fmt(trail[0][0])}, {K.fmt(trail[0][1])})"
          f"   <- trail input")
    print(f"(dX^{R_IN+R}, dY^{R_IN+R}) = ({K.fmt(a)}, {K.fmt(b)})   <- (a, b), forced by "
          f"the ciphertext pairing")

    print("\n===== key-free filters =====")
    f0, f1 = wt(e["m0"]), wt(e["m1"])
    print(f"dX^0 & {K.fmt(e['m0'])} = {K.fmt(e['v0'])}   wt = {f0}")
    print(f"dX^1 & {K.fmt(e['m1'])} = {K.fmt(e['v1'])}   wt = {f1}")
    print(f"F = {e['F']}")

    print("\n===== round-key masks =====")
    for t in sorted(e["masks"]):
        m = e["masks"][t]
        print(f"M_K{t-1} (sk^{t-1}) = {K.fmt(m)}   wt = {wt(m)}")
    k = e["k"]
    print(f"k = {k}")

    # ------------------------------------------------------------ complexities
    N = N_PAIRS
    resid3 = RESID3
    resid2 = BLOCK - f0 - f1 - resid3

    # 段階的な早期打ち切り: 第 1 段階で sk^0 の M_K0^(1) を推測して ΔX^2 の
    # resid2 ビットを検査し, 第 2 段階で残りの sk^0 と sk^1 を推測する。
    sA = (N - f0 - f1) + STAGE1_BITS
    survA = sA - resid2
    sB = survA + (wt(e["masks"][1]) - STAGE1_BITS) + wt(e["masks"][2])

    # 時間はすべて 23 段暗号化の 2/23 を単位に数える
    u = math.log2(2 / total_rounds)
    t_data = u + (N + 1)                    # Steps 1, 2 : データ生成
    t_dx1 = u + (N - f0)                    # Step 3     : ΔX^1 の再計算
    t_st1 = u + sA                          # Stage 1
    t_st2 = u + sB                          # Stage 2
    t_fin = max(LOG_E_WRONG, 0.0) + (KEYLEN - k)    # Step 5 : 残る鍵の全数探索
    T = math.log2(2 ** t_data + 2 ** t_dx1 + 2 ** t_st1 + 2 ** t_st2 + 2 ** t_fin)
    T_kr = math.log2(2 ** t_dx1 + 2 ** t_st1 + 2 ** t_st2 + 2 ** t_fin)
    mem = math.log2(2 * 2 ** (N - f0 - f1) * 16)    # 対あたり 2 テキスト x 16 バイト

    lam_r = 2.0 ** (N - W)
    lam_w = 2.0 ** (N - BLOCK)
    p_succ = poisson_ge(lam_r, TAU)
    log_e_wrong = (-lam_w / math.log(2) + TAU * math.log2(lam_w)
                   - math.log2(math.factorial(TAU))) + k

    print(f"\n===== complexities ({total_rounds}-round attack) =====")
    print(f"pairs N                = 2^{N}")
    print(f"data D                 = 2^{N + 1} chosen ciphertexts")
    print(f"right pairs            = 2^{N - W} = {2 ** (N - W)}")
    print(f"after Eq.(11)          = 2^{N - f0}")
    print(f"after Eq.(12)          = 2^{N - f0 - f1}")
    print(f"wrong-class expectation= 2^{N - BLOCK}")
    print(f"residual conditions    = {f0} + {f1} + {resid2} + {resid3} = "
          f"{f0 + f1 + resid2 + resid3}")
    print(f"stage 1                : 2^{sA} -> 2^{survA} survivors")
    print(f"stage 2                : 2^{sB} -> 2^{sB - resid3} counter increments")
    print(f"time  Steps 1,2 = 2^{t_data:.2f} | Step 3 = 2^{t_dx1:.2f} | "
          f"Stage 1 = 2^{t_st1:.2f} | Stage 2 = 2^{t_st2:.2f} | Step 5 = 2^{t_fin:.2f}")
    print(f"time  T                = 2^{T:.2f}")
    print(f"time  T (key recovery) = 2^{T_kr:.2f}")
    print(f"memory M               = 2^{mem:.2f} bytes")
    print(f"threshold tau          = {TAU}")
    print(f"  Pr[C_right >= tau]   = {p_succ:.5f}")
    print(f"  E[# wrong >= tau]    = 2^{log_e_wrong:.1f}")

    # -------------------------------------------------- agreement with the paper
    print("\n===== agreement with paper/SAND_DC_IEICE (Table 13, Appendix) =====")
    check("the selected characteristic is the one of Table 12",
          (trail[0] == PAPER_TRAIL_IN and trail[R] == PAPER_TRAIL_OUT))
    check(f"wt(mask of Eq.(11)) = {PAPER['f0']}", f0 == PAPER["f0"], str(f0))
    check(f"wt(mask of Eq.(12)) = {PAPER['f1']}", f1 == PAPER["f1"], str(f1))
    check(f"wt(M_K0) = {PAPER['M_K0']}", wt(e["masks"][1]) == PAPER["M_K0"],
          str(wt(e["masks"][1])))
    check(f"wt(M_K1) = {PAPER['M_K1']}", wt(e["masks"][2]) == PAPER["M_K1"],
          str(wt(e["masks"][2])))
    check(f"k = {PAPER['k']}", k == PAPER["k"], str(k))
    check("residual conditions sum to the 128-bit condition",
          f0 + f1 + resid2 + resid3 == BLOCK)
    check(f"staged counting reproduces 2^{k} x 2^{N - BLOCK} increments",
          sB - resid3 == k + (N - BLOCK), f"2^{sB - resid3} vs 2^{k + N - BLOCK}")
    check(f"data = 2^{PAPER['data']}", N + 1 == PAPER["data"], f"2^{N + 1}")
    check(f"T = 2^{PAPER['T']}", abs(T - PAPER["T"]) < 0.1, f"2^{T:.2f}")
    check(f"T (key recovery) = 2^{PAPER['T_kr']}", abs(T_kr - PAPER["T_kr"]) < 0.1,
          f"2^{T_kr:.2f}")
    check(f"memory = 2^{PAPER['mem']} bytes", abs(mem - PAPER["mem"]) < 0.1,
          f"2^{mem:.2f}")
    check(f"success probability = {PAPER['p_succ']}",
          abs(p_succ - PAPER["p_succ"]) < 1e-4, f"{p_succ:.5f}")
    check("no wrong class is expected to survive the threshold", log_e_wrong < 0,
          f"2^{log_e_wrong:.1f}")
    check(f"data 2^{N + 1} < full codebook 2^{BLOCK}", N + 1 < BLOCK)
    check(f"time 2^{T:.2f} < exhaustive search 2^{KEYLEN}", T < KEYLEN)

    print("\n" + "=" * 78)
    print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
    print("=" * 78)
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
