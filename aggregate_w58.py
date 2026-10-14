#!/usr/bin/env python3
"""
列挙結果 (trails_parts/ と trails_w58_fine/) を集約・重複除去し, 統計を出したうえで
SAT モデルとは独立な転送行列計算で全件を検証する。

出力
  * 重み 58 の 12 段特性の総数
  * 重みベクトルごとの内訳
  * 相異なる入力差分 / 出力差分 / (入力, 出力) 対の個数
  * 特性数が多い (入力, 出力) 対 = 差分クラスタリングが強い差分
  * 各特性の全ラウンド遷移が転送行列計算で実現可能か
"""
import glob
import importlib.util
import math
import multiprocessing as mp
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

spec = importlib.util.spec_from_file_location(
    "transfer_matrix", os.path.join(HERE, "transfer_matrix.py"))
V = importlib.util.module_from_spec(spec)
sys.modules["transfer_matrix"] = V
spec.loader.exec_module(V)

NUM_ROUNDS = 12
DIRS = [os.path.join(HERE, "trails_parts"), os.path.join(HERE, "trails_w58_fine")]
OUT = os.path.join(HERE, "trails_w58_all.txt")


def parse(path):
    trails, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# characteristic"):
                if cur:
                    trails.append(cur)
                cur = []
            elif line and not line.startswith("#"):
                p = line.split()
                if len(p) == 5:
                    cur.append((int(p[1] + p[2], 16), int(p[3] + p[4], 16)))
    if cur:
        trails.append(cur)
    return [t for t in trails if len(t) == NUM_ROUNDS + 1]


def load_all():
    seen, out = set(), []
    for d in DIRS:
        for path in sorted(glob.glob(os.path.join(d, "*.txt"))):
            for t in parse(path):
                k = tuple(t)
                if k not in seen:
                    seen.add(k)
                    out.append(t)
    return out


def check_trail(trail):
    """全ラウンド遷移の厳密重み (32 - log2 count). 不可能な遷移があれば None"""
    w = []
    for r in range(len(trail) - 1):
        c = V.round_transition_count(*trail[r], *trail[r + 1])
        if c == 0:
            return None
        w.append(32.0 - math.log2(c))
    return w


def fmt(word):
    return f"{word >> 16:04X} {word & 0xFFFF:04X}"


def main():
    trails = load_all()
    n = len(trails)
    print(f"SAND-64 : 12-round differential characteristics with p = 2^-58")
    print(f"  (weight 57 is UNSAT, so 58 is the minimum and every trail below is optimal)\n")
    print(f"total distinct characteristics : {n}")

    ins = Counter((t[0][0], t[0][1]) for t in trails)
    outs = Counter((t[-1][0], t[-1][1]) for t in trails)
    pairs = Counter(((t[0][0], t[0][1]), (t[-1][0], t[-1][1])) for t in trails)
    print(f"distinct input  differences    : {len(ins)}")
    print(f"distinct output differences    : {len(outs)}")
    print(f"distinct (input, output) pairs : {len(pairs)}")

    # 中間状態の活性ニブル数の分布など, 構造の概観
    print("\n-- number of trails per (input,output) pair --")
    hist = Counter(pairs.values())
    for k in sorted(hist):
        print(f"   {k} trail(s) : {hist[k]} pairs")

    print("\n-- (input, output) pairs carrying the most weight-58 trails (top 15) --")
    for (i, o), c in pairs.most_common(15):
        print(f"   dP=({fmt(i[0])},{fmt(i[1])}) -> dC=({fmt(o[0])},{fmt(o[1])}) : {c} trails"
              f"  => differential p >= 2^-{58 - math.log2(c):.2f}")

    print("\n-- independent transfer-matrix verification of every characteristic --", flush=True)
    with mp.Pool(processes=15) as pool:
        res = pool.map(check_trail, trails, chunksize=8)
    bad = [i for i, r in enumerate(res) if r is None]
    print(f"   characteristics with an impossible round transition : {len(bad)}")
    sums = Counter(round(sum(r), 6) for r in res if r is not None)
    print("   exact per-round (intra-round clustered) weight sum:")
    for k in sorted(sums):
        print(f"     2^-{k:g} : {sums[k]} characteristics")

    with open(OUT, "w") as f:
        f.write(f"# SAND-64 : all {n} 12-round differential characteristics with p = 2^-58\n")
        f.write("# 58 is the minimum weight over all 12-round characteristics (weight <= 57 is UNSAT).\n")
        f.write("# Each block lists r, dX^r, dY^r for r = 0..12 (dX^0||dY^0 = plaintext difference).\n")
        for idx, t in enumerate(trails, 1):
            f.write(f"\n# characteristic {idx}\n")
            for r, (dx, dy) in enumerate(t):
                f.write(f"{r:>3} {fmt(dx)} {fmt(dy)}\n")
    print(f"\nall {n} characteristics written to {OUT}")


if __name__ == "__main__":
    main()
