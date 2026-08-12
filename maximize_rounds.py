#!/usr/bin/env python3
"""
SAND-128 (ブロック 128 bit, 鍵 128 bit) で攻撃可能段数を最大化する。

攻撃段数 = r_in + R + r_out
  R      : 識別子の段数。D ~ 2^(W+5) < 2^128 より W <= 123。B[19]=120, B[20]>=126
           なので R = 19 が上限。
  r_out  : 論文と同じく暗号文側 1 段 (最終ラウンド鍵に依存しない対応付けで無料)
  r_in   : 平文側に付けられる段数。鍵非依存フィルタ F(r_in) が効く範囲で伸ばせる。

成立条件 (対数 N = 2^(W+4), 正対 2^4 個):
    誤対の残存 / 鍵クラス数 = 2^(W+4-F-k) < 2^4   <=>   F + k > W

F を大きく保つには平文側境界を疎にする必要があるので, R=19 で境界の活性ニブル数を
最小化する特性を探す。

    python3 maximize_rounds.py [workers]
"""
import multiprocessing as mp
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import keyrecovery_sand128 as K
import filters_sand128 as F
import matsui_sand128 as M

ROUNDS = 19
BLOCK = 128
KEY = 128
MARGIN = 4          # 対数 N = 2^(W+MARGIN) -> 正対 2^MARGIN 個
R_OUT = 1


def job(args):
    weight, n_pt, timeout = args
    bounds = M.load_bounds()
    m, variables = K.build(ROUNDS, weight, bounds, max_nib_pt=n_pt)
    sat, trail, rw, dt = K.solve(m, variables, ROUNDS, f"mx_{weight}_{n_pt}", timeout)
    if not sat:
        return weight, n_pt, None, dt
    return weight, n_pt, trail, dt


def evaluate(trail, weight, max_r_in=8):
    """境界の疎さからフィルタ減衰を計算し, 成立する最大の r_in を求める"""
    # 平文側境界 (挿入後の第 s 段) = (ΔY^R, ΔX^R)
    dx_s = trail[ROUNDS][1]        # ΔY^R
    dy_s = trail[ROUNDS][0]        # ΔX^R
    seq = F.truncated_backward(dx_s, dy_s, max_r_in, 16, K.P_BOX)

    rows = []
    for r_in in range(2, max_r_in + 1):
        dx0, _ = seq[r_in]
        dx1, _ = seq[r_in - 1]
        k0 = sum(1 for b in dx0 if b is not None)
        k1 = sum(1 for b in dx1 if b is not None)
        f_bits = k0 + k1
        # 鍵推測なし (k=0) で成立するか / 必要な k
        need_k = max(0, weight - f_bits + 1)
        rows.append((r_in, k0, k1, f_bits, need_k))
    return dx_s, dy_s, rows


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 5400.0

    # D = 2^(W+5) < 2^128 -> W <= 123
    weights = [120, 122]
    n_pts = [1, 2, 3, 4, 5, 6]
    jobs = [(w, n, timeout) for w in weights for n in n_pts]
    print(f"R={ROUNDS}, block {BLOCK}, key {KEY}; searching sparsest plaintext-side boundary")
    print(f"{len(jobs)} jobs, {workers} workers, timeout {timeout:.0f}s", flush=True)

    t0 = time.time()
    with mp.Pool(processes=workers) as pool:
        results = pool.map(job, jobs, chunksize=1)

    print(f"\nsearch done ({time.time()-t0:.0f}s)")
    print(f"{'weight':>7} {'n_pt':>5} {'result':>10} {'time':>8}")
    found = {}
    for weight, n_pt, trail, dt in results:
        tag = "SAT" if trail else "UNSAT"
        print(f"{weight:>7} {n_pt:>5} {tag:>10} {dt:>7.0f}s")
        if trail and (weight not in found or n_pt < found[weight][0]):
            found[weight] = (n_pt, trail)

    best = None
    for weight in sorted(found):
        n_pt, trail = found[weight]
        exact_w, total = K.verify_trail(trail, f"R19 w{weight}")
        dx_s, dy_s, rows = evaluate(trail, weight)
        print(f"\n===== weight {weight} (transfer-matrix verified: {total}), "
              f"plaintext-side {n_pt} nibbles =====")
        print(f"  plaintext-side boundary: dX^s = {K.fmt(dx_s)}, dY^s = {K.fmt(dy_s)}")
        print(f"  {'r_in':>5} {'|dX0|':>6} {'|dX1|':>6} {'F':>5} {'need k':>7} {'rounds':>7}")
        for r_in, k0, k1, f_bits, need_k in rows:
            rounds = r_in + ROUNDS + R_OUT
            print(f"  {r_in:>5} {k0:>6} {k1:>6} {f_bits:>5} {need_k:>7} {rounds:>7}")
            cand = (rounds, need_k, weight, r_in, trail)
            if need_k <= 0 and (best is None or rounds > best[0]):
                best = cand

    if best:
        rounds, need_k, weight, r_in, trail = best
        print(f"\n### key-free filter alone suffices: {rounds} rounds "
              f"(r_in={r_in} + R={ROUNDS} + r_out={R_OUT}), weight {weight}")
        print(f"    D = 2^{weight + MARGIN + 1} chosen ciphertexts, "
              f"N = 2^{weight + MARGIN} pairs, right pairs = 2^{MARGIN}")
    return found


if __name__ == "__main__":
    main()
