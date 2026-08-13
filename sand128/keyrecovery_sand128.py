#!/usr/bin/env python3
"""
SAND-128 : 選択暗号文設定の鍵回復攻撃に向けた差分特性の探索と検証。

CANDAR W 2026 の SAND-64 論文 (Differential Attack on AND-RX-based Lightweight
Block Cipher SAND-64) の構成を SAND-128 に移す。

1. 暗号文側からの挿入 (選択暗号文設定)
   R 段の識別子を攻撃の第 s 段から第 s+R 段に置く。復号方向のラウンド関数は
   暗号化のラウンド関数の左右半分を入れ替えたものに等しいので, 順方向の特性
   (ΔX^t, ΔY^t) は攻撃側の番号付けでは

       (ΔX^{s+R-t}, ΔY^{s+R-t}) = (ΔY^t, ΔX^t)

   となる。実際, U^t = ΔY^{R-t}, V^t = ΔX^{R-t} と置くと r = R-t-1 に対して
       V^{t+1} = ΔX^{r}   = U^t
       U^{t+1} = ΔY^{r}   = ΔX^{r+1} ⊕ P(ΔG0 ⊕ ΔG1) = V^t ⊕ P(F(U^t))
   であり順方向と同一の関係式を満たす。ラウンド重みは順序が反転する。
   論文の (ΔX^15, ΔY^15) = (00182C28, 00404000) が Table III の
   (ΔX^0, ΔY^0) = (00404000, 00182C28) の swap であることと一致する。

2. 手順 4 : 鍵回復コストを目的関数に入れる
   鍵推測に必要なビット数を決めるのは差分の「活性ニブル (SSb 位置) の個数」で
   あって, 活性ビット数ではない。1 活性ニブルが 1 サブ鍵ニブルの推測に対応する
   ためである。各ニブル j に指示変数 a_j = OR(4ビット) を導入し, Σ a_j に
   カーディナリティ制約を課して端の疎さを最適化する。

3. 転送行列による厳密検証 (SAT モデルに依存しない独立検証)
   ラウンド遷移を実現する 64 ビットラウンド入力の個数を, 16 個の SSb 位置に
   わたる転送行列積で厳密に数える。G0 と G1 が同じニブルを共有する制約
   (joint DDT) を全位置で強制するので, 独立性の仮定を一切使わない。

使い方
    python3 keyrecovery_sand128.py trail <rounds> <weight> [max_nib_pt] [max_nib_ct]
    python3 keyrecovery_sand128.py pareto <rounds> [workers]
"""
import importlib.util
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matsui_sand128 as M

sand128 = M.sand128
PER_ROUND = M.PER_ROUND
P_BOX = [14, 15, 8, 9, 2, 3, 12, 13, 6, 7, 0, 1, 10, 11, 4, 5]
TMP_DIR = "/dev/shm/sand128_kr"


# ==========================================================================
# ビット/ニブルの取り扱い  (index 0 = MSB of the 64-bit half)
# ==========================================================================
def bit(word, i):
    return (word >> (63 - i)) & 1


def nibble(word, j):
    """SSb ニブル j = ビット (j, j+16, j+32, j+48), 先頭が MSB"""
    return (bit(word, j) << 3) | (bit(word, j + 16) << 2) | (bit(word, j + 32) << 1) | bit(word, j + 48)


def active_nibbles(word):
    return [j for j in range(16) if nibble(word, j)]


def apply_p(word):
    """out index i+16k <- in index P_BOX[i]+16k"""
    out = 0
    for k in range(4):
        for i in range(16):
            out |= bit(word, P_BOX[i] + 16 * k) << (63 - (i + 16 * k))
    return out


def apply_p_inv(word):
    out = 0
    for k in range(4):
        for i in range(16):
            out |= bit(word, i + 16 * k) << (63 - (P_BOX[i] + 16 * k))
    return out


# ==========================================================================
# 転送行列による厳密なラウンド遷移個数
# ==========================================================================
_DDT = None


def ddt():
    global _DDT
    if _DDT is None:
        _DDT = sand128.SAND_128_SAT_Model.ssb_joint_ddt()
    return _DDT


def round_transition_count(dx_prev, dy_prev, dx_cur):
    """
    ΔX^{r-1} = dx_prev, ΔY^{r-1} = dy_prev, ΔX^r = dx_cur を実現する
    64 ビットラウンド入力の個数を厳密に数える。

    ΔX^r = ΔY^{r-1} ⊕ P(ΔG0 ⊕ ΔG1) より T = P^{-1}(ΔX^r ⊕ ΔY^{r-1}) = ΔG0 ⊕ ΔG1。
    入力ニブル m に対し (d0_m, d1_m) = (G0 の出力差分ニブル m, G1 の出力差分ニブル m-1)
    で, DDT_SSb[din_m][(d0_m<<4)|d1_m] がその個数を与える。
    ΔG1 のニブル j は入力ニブル j+1 由来なので, 位置 j の等式は
        T_j = d0_j ⊕ d1_{j+1}
    となり, d0 を 1 つ先へ送る巡回的な DP (転送行列) になる。
    """
    T = apply_p_inv(dx_cur ^ dy_prev)
    tab = ddt()
    din = [nibble(dx_prev, m) for m in range(16)]
    tj = [nibble(T, j) for j in range(16)]

    total = 0
    # d0_0 の値ごとに鎖を開き, 最後に巡回条件で閉じる
    for d0_first in range(16):
        # state: d0_m の値 -> 場合の数
        state = {d0_first: 0}
        # m = 0 の (d0_0, d1_0) を選ぶ。d1_0 は位置 15 の等式で使うので後回しにできない
        # -> d1_0 も固定して鎖を回す
        for d1_first in range(16):
            c0 = tab[din[0]][(d0_first << 4) | d1_first]
            if c0 == 0:
                continue
            cur = {d0_first: c0}
            for m in range(1, 16):
                nxt = {}
                for d0_prev, ways in cur.items():
                    need_d1 = tj[m - 1] ^ d0_prev        # T_{m-1} = d0_{m-1} ⊕ d1_m
                    row = tab[din[m]]
                    for d0 in range(16):
                        c = row[(d0 << 4) | need_d1]
                        if c:
                            nxt[d0] = nxt.get(d0, 0) + ways * c
                cur = nxt
                if not cur:
                    break
            # 巡回条件: T_15 = d0_15 ⊕ d1_0  (d1_0 = d1_first)
            for d0_last, ways in cur.items():
                if tj[15] ^ d0_last == d1_first:
                    total += ways
    return total


def verify_trail(trail, label=""):
    """
    trail = [(dx, dy), ...] (r = 0..R) の各ラウンド遷移を転送行列で厳密検証する。
    返り値: (ラウンド重みのリスト, 総重み) / 不可能遷移があれば例外
    """
    weights = []
    for r in range(1, len(trail)):
        dx_prev, dy_prev = trail[r - 1]
        dx_cur, dy_cur = trail[r]
        if dy_cur != dx_prev:
            raise AssertionError(f"{label} round {r}: ΔY^r != ΔX^(r-1)")
        cnt = round_transition_count(dx_prev, dy_prev, dx_cur)
        if cnt == 0:
            raise AssertionError(f"{label} round {r}: impossible transition (count 0)")
        if cnt & (cnt - 1):
            raise AssertionError(f"{label} round {r}: count {cnt} is not a power of two")
        w = 64 - cnt.bit_length() + 1
        weights.append(w)
    return weights, sum(weights)


# ==========================================================================
# SAT 探索 : 総重み <= W, Matsui 枝刈り, 端の活性ニブル数の制約
# ==========================================================================
def add_nibble_indicators(m, variables, r, half):
    """半分 (half='X' -> ΔX^r, 'Y' -> ΔY^r) のニブル活性指示変数 16 個を作る"""
    base = 0 if half == "X" else 64
    ind = []
    for j in range(16):
        m.aux_var_counter += 1
        a = m.aux_var_counter
        bits = [variables[f"x_{base + idx}_{r}"] for idx in (j, j + 16, j + 32, j + 48)]
        for b in bits:
            m.clauses.append([-b, a])          # bit -> active
        m.clauses.append([-a] + bits)          # active -> some bit
        ind.append(a)
    return ind


def build(rounds, weight, bounds, max_nib_pt=None, max_nib_ct=None):
    """
    rounds 段・総重み <= weight のモデル。

    max_nib_pt : 攻撃の平文側境界 (= 特性の ΔX^R と ΔX^{R-1}) の活性ニブル数の上限
    max_nib_ct : 攻撃の暗号文側境界 (= 特性の ΔX^0 と ΔY^0) の活性ニブル数の上限
    """
    m = M.CounterModel()
    variables = m.generate_variables(rounds)
    m.Search_Differential_Distingusher(rounds, weight, variables)
    m.clauses.append([variables[f"x_{i}_0"] for i in range(128)])

    # --- Matsui 型の下界枝刈り (counter_s を上書きする前に入れる) ---
    for k in range(1, rounds):
        b = bounds.get(rounds - k)
        if b is None:
            continue
        limit = weight - b
        if limit < 0:
            m.clauses.append([])
            return m, variables
        if (PER_ROUND * k, limit + 1) in m.counter_s:
            m.clauses.append([-m.counter_s[(PER_ROUND * k, limit + 1)]])

    # --- 端の活性ニブル数の制約 ---
    if max_nib_ct is not None:
        ind = add_nibble_indicators(m, variables, 0, "X") + \
              add_nibble_indicators(m, variables, 0, "Y")
        m.vars = ind
        m.sequential_encoding_method(max_nib_ct)
    if max_nib_pt is not None:
        ind = add_nibble_indicators(m, variables, rounds, "X") + \
              add_nibble_indicators(m, variables, rounds - 1, "X")
        m.vars = ind
        m.sequential_encoding_method(max_nib_pt)

    return m, variables


def solve(m, variables, rounds, tag, timeout):
    if any(len(c) == 0 for c in m.clauses):
        return False, None, None, 0.0
    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"{tag}.cnf")
    with open(cnf, "w") as f:
        f.write(f"p cnf {m.aux_var_counter} {len(m.clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in m.clauses))
    t0 = time.time()
    sat, pos = M.run_kissat(m.kissat_path, cnf, timeout)
    dt = time.time() - t0
    if os.path.exists(cnf):
        os.remove(cnf)
    if not sat:
        return sat, None, None, dt

    trail = []
    for r in range(rounds + 1):
        dx = 0
        dy = 0
        for i in range(64):
            if variables[f"x_{i}_{r}"] in pos:
                dx |= 1 << (63 - i)
            if variables[f"x_{i + 64}_{r}"] in pos:
                dy |= 1 << (63 - i)
        trail.append((dx, dy))
    rw = [sum(1 for i in range(PER_ROUND) if variables[f"ssb_w_{i}_{r}"] in pos)
          for r in range(1, rounds + 1)]
    return True, trail, rw, dt


# ==========================================================================
# 暗号文側からの挿入
# ==========================================================================
def insert_backward(trail, s):
    """
    順方向の特性 trail (t = 0..R) を, 攻撃の第 s 段から第 s+R 段に
    暗号文側から挿入した形に変換する。
        (ΔX^{s+R-t}, ΔY^{s+R-t}) = (ΔY^t, ΔX^t)
    返り値: {攻撃ラウンド番号: (ΔX, ΔY)}
    """
    R = len(trail) - 1
    return {s + R - t: (dy, dx) for t, (dx, dy) in enumerate(trail)}


def fmt(word):
    return f"{word >> 32:08X}{word & 0xFFFFFFFF:08X}"


def report(trail, rounds, rw, s, label):
    exact_w, exact_total = verify_trail(trail, label)
    print(f"\n=== {label} : forward characteristic (R={rounds}) ===")
    print(f"{'t':>3} {'dX^t':>17} {'dY^t':>17} {'w':>3} {'exact':>6} {'log2(p)':>8}  act(dX)")
    acc = 0
    for t, (dx, dy) in enumerate(trail):
        if t == 0:
            print(f"{t:>3} {fmt(dx):>17} {fmt(dy):>17} {'':>3} {'':>6} {0:>8}  {active_nibbles(dx)}")
        else:
            acc += exact_w[t - 1]
            flag = "" if rw[t - 1] == exact_w[t - 1] else f" (SAT said {rw[t-1]})"
            print(f"{t:>3} {fmt(dx):>17} {fmt(dy):>17} {rw[t-1]:>3} {exact_w[t-1]:>6} {-acc:>8}"
                  f"  {active_nibbles(dx)}{flag}")
    print(f"SAT model total weight   : {sum(rw)}")
    print(f"transfer-matrix total    : {exact_total}   ->  p = 2^-{exact_total}")
    if sum(rw) != exact_total:
        print("!! SAT model and transfer matrix disagree")

    ins = insert_backward(trail, s)
    R = rounds
    print(f"\n=== {label} : inserted from the ciphertext side, rounds {s}..{s+R} ===")
    print("(chosen-ciphertext setting; halves are swapped, as in the CANDAR paper)")
    print(f"{'round':>5} {'dX^round':>17} {'dY^round':>17}  act(dX)  act(dY)")
    for a in range(s, s + R + 1):
        dx, dy = ins[a]
        print(f"{a:>5} {fmt(dx):>17} {fmt(dy):>17}  {len(active_nibbles(dx)):>5}  "
              f"{len(active_nibbles(dy)):>6}")

    ct = ins[s + R]
    pt = ins[s]
    print(f"\nciphertext-side boundary (round {s+R}): "
          f"a = dX = {fmt(ct[0])} ({len(active_nibbles(ct[0]))} nibbles), "
          f"b = dY = {fmt(ct[1])} ({len(active_nibbles(ct[1]))} nibbles)")
    print(f"plaintext-side boundary  (round {s}): "
          f"dX = {fmt(pt[0])} ({len(active_nibbles(pt[0]))} nibbles), "
          f"dY = {fmt(pt[1])} ({len(active_nibbles(pt[1]))} nibbles)")
    n_pt = len(active_nibbles(pt[0])) + len(active_nibbles(pt[1]))
    n_ct = len(active_nibbles(ct[0])) + len(active_nibbles(ct[1]))
    print(f"key-recovery proxy: plaintext-side active nibbles = {n_pt}, "
          f"ciphertext-side = {n_ct}")

    kc = key_cost(trail, s)
    print(f"\n--- required round-key bits (paper's M_K0 / M_K1 analogue) ---")
    print(f"  test dX^{s}   = {fmt(kc['dX_s'])}  ({kc['act_dX_s']} active nibbles)"
          f"  -> mask {fmt(kc['mask_sk_s_minus_1'])}  wt = {kc['wt1']}")
    print(f"  test dY^{s}   = {fmt(kc['dY_s'])}  ({kc['act_dY_s']} active nibbles)"
          f"  -> mask {fmt(kc['mask_sk_s_minus_2'])}  wt = {kc['wt0']}")
    print(f"  total signal-bearing round-key bits = {kc['wt1'] + kc['wt0']}"
          f"   -> at most 2^{kc['wt1'] + kc['wt0']} partial-key classes")
    print(f"  per-nibble derivative支持: "
          + ", ".join(f"n{j}:d=0x{nibble(kc['dX_s'], j):X}->{len(derivative_support(nibble(kc['dX_s'], j)))}b"
                      for j in active_nibbles(kc['dX_s'])))
    return exact_total, n_pt, n_ct, kc


# ==========================================================================
# 論文 Appendix と同じサブ鍵マスクの算出
# ==========================================================================
def ssb_value(u):
    G = sand128.SAND_128_SAT_Model
    return (G.g0_value(u) << 4) | G.g1_value(u)


def derivative_support(delta):
    """
    入力差分 delta のニブルに対し, 導関数 D_delta SSb(u) = SSb(u) ^ SSb(u ^ delta)
    が実際に依存する入力ビットの集合 (3 = msb, 0 = lsb) を返す。

    論文 Appendix A と同じ計算。delta = 0x4 のとき
        D_delta N0(u) = (u3u1, 1, 0, u3),  D_delta N1(u) = (0, 1, u0, 0)
    なので依存先は {u3, u1, u0} の 3 ビットで, u2 には依存しない。
    """
    if delta == 0:
        return set()

    def D(u):
        return ssb_value(u) ^ ssb_value(u ^ delta)

    dep = set()
    for b in range(4):
        if any(D(u) != D(u ^ (1 << b)) for u in range(16)):
            dep.add(b)
    return dep


def required_mask(diff_word, ncols=16):
    """
    差分 diff_word を検査するために計算が必要な状態座標のマスクを返す。

    語は (u3 行)||(u2 行)||(u1 行)||(u0 行) の順に並ぶ (論文の x{3}||x{2}||x{1}||x{0})。
    行 b の列 j はビット位置 (3-b)*ncols + j にある。
    """
    nbits = 4 * ncols
    mask = 0
    for j in range(ncols):
        if ncols == 16:
            d = nibble(diff_word, j)
        else:
            d = 0
            for k, b in enumerate((3, 2, 1, 0)):
                if (diff_word >> (nbits - 1 - (k * ncols + j))) & 1:
                    d |= 1 << b
        for b in derivative_support(d):
            mask |= 1 << (nbits - 1 - ((3 - b) * ncols + j))
    return mask


def validate_against_paper():
    """
    CANDAR 論文 (SAND-64) の Appendix A の数値でマスク算出を検証する。
    ΔX^2 = 008C0000 に対し M_K1 = 8C008C8C, wt = 9 が得られるはず。
    G0, G1 は SAND-64 と SAND-128 で同一なのでこの検証はそのまま通用する。
    """
    dep = derivative_support(0x4)
    assert dep == {3, 1, 0}, dep
    m = required_mask(0x008C0000, ncols=8)
    assert m == 0x8C008C8C, f"{m:08X}"
    assert bin(m).count("1") == 9
    print("paper check: derivative_support(0x4) = {u3,u1,u0}  (3 bits)")
    print(f"paper check: required_mask(dX2=008C0000) = {m:08X}, wt = {bin(m).count('1')}"
          f"   -> matches M_K1 = 8C008C8C, wt 9")
    return True


def key_cost(trail, s):
    """
    暗号文側から挿入したときの平文側境界について, 必要な丸鍵ビット数を見積る。

    論文と同じ構成 (平文側に 3 段, 暗号文側に 1 段) を想定する。識別子の先頭が
    攻撃の第 s 段なので, 検査すべき差分は攻撃側の ΔX^s と ΔX^{s-1} である。
    """
    ins = insert_backward(trail, s)
    dx_s, dy_s = ins[s]

    # 攻撃の第 s 段の差分を検査するために必要な X^{s-1} の座標
    m1 = required_mask(dx_s)
    # ΔY^s = ΔX^{s-1} なので, 1 段手前の検査に必要な座標
    m0 = required_mask(dy_s)
    return {
        "dX_s": dx_s, "dY_s": dy_s,
        "mask_sk_s_minus_1": m1, "wt1": bin(m1).count("1"),
        "mask_sk_s_minus_2": m0, "wt0": bin(m0).count("1"),
        "act_dX_s": len(active_nibbles(dx_s)),
        "act_dY_s": len(active_nibbles(dy_s)),
    }


# ==========================================================================
# Pareto 走査
# ==========================================================================
def _job(args):
    rounds, weight, n_pt, bounds, timeout = args
    m, variables = build(rounds, weight, bounds, max_nib_pt=n_pt)
    sat, trail, rw, dt = solve(m, variables, rounds, f"p_{rounds}_{weight}_{n_pt}", timeout)
    return weight, n_pt, sat, dt, (trail, rw) if sat else None


def pareto(rounds, workers, timeout=1800.0):
    bounds = M.load_bounds()
    base = bounds.get(rounds)
    if base is None:
        print(f"B[{rounds}] が bounds_sand128.json にありません")
        return
    print(f"pareto scan: R={rounds}, B[{rounds}]={base}, workers={workers}")
    jobs = []
    for weight in range(base, base + 13, 2):
        for n_pt in range(2, 15):
            jobs.append((rounds, weight, n_pt, bounds, timeout))
    with mp.Pool(processes=workers) as pool:
        results = pool.map(_job, jobs, chunksize=1)

    grid = {}
    best = {}
    for weight, n_pt, sat, dt, payload in results:
        grid[(weight, n_pt)] = sat
        if sat and (weight not in best or n_pt < best[weight][0]):
            best[weight] = (n_pt, payload)

    print(f"\nSAT (o) / UNSAT (.) / timeout (?)   rows: total weight, cols: plaintext-side nibbles")
    header = "      " + "".join(f"{n:>3}" for n in range(2, 15))
    print(header)
    for weight in range(base, base + 13, 2):
        row = "".join("  o" if grid.get((weight, n)) else ("  ?" if grid.get((weight, n)) is None else "  .")
                      for n in range(2, 15))
        print(f"w={weight:>3} {row}")

    print("\nminimum plaintext-side active nibbles per weight budget:")
    for weight in sorted(best):
        n_pt, payload = best[weight]
        print(f"  weight {weight}: {n_pt} nibbles")
    return best


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "trail"
    bounds = M.load_bounds()

    if cmd == "trail":
        rounds = int(sys.argv[2])
        weight = int(sys.argv[3])
        n_pt = int(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != "-" else None
        n_ct = int(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != "-" else None
        s = int(os.environ.get("SAND128_S", "3"))     # 識別子を置く先頭ラウンド
        m, variables = build(rounds, weight, bounds, n_pt, n_ct)
        print(f"R={rounds}, weight<={weight}, max_nib_pt={n_pt}, max_nib_ct={n_ct}")
        print(f"CNF: {m.aux_var_counter} vars, {len(m.clauses)} clauses", flush=True)
        sat, trail, rw, dt = solve(m, variables, rounds, f"t_{rounds}_{weight}", 7200.0)
        if not sat:
            print(f"UNSATISFIABLE ({dt:.1f}s)")
            return
        print(f"SATISFIABLE ({dt:.1f}s)")
        report(trail, rounds, rw, s, f"R={rounds} w<={weight}")

    elif cmd == "pareto":
        rounds = int(sys.argv[2])
        workers = int(sys.argv[3]) if len(sys.argv) > 3 else 14
        pareto(rounds, workers)


if __name__ == "__main__":
    main()
