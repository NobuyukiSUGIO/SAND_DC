#!/usr/bin/env python3
"""
enum_by_wvec.py の phase B をさらに細かく分割して並列度を上げる。

分割の第 2 軸: 「ラウンド 1 で活性な SSb の位置集合 S」
重みベクトルにより ラウンド 1 の重みは 2 (SSb 1 個が重み 2) または
4 (SSb 2 個が重み 2; SSb の重みは 0,2,3 のいずれかなので 4 = 2+2 しかない)
に固定されているので, S の大きさは 1 または 2 である。
S は互いに素で全体を覆うから, (重みベクトル, S) が正しい分割になる。

既に trails_w58/ に見つかっている特性はブロッキング節として読み込むので,
これまでの計算結果は無駄にならない。
"""
import glob
import itertools
import json
import multiprocessing as mp
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from enumerate_w58 import (key_vars, write_cnf, run_kissat, extract_trail, fmt)
import enum_by_wvec as W

NUM_ROUNDS = W.NUM_ROUNDS
OUT_DIR = os.path.join(HERE, "trails_w58_fine")
PARTS_DIR = os.path.join(HERE, "trails_parts")
# CNF は毎回書き直すので tmpfs (RAM) に置く。ディスク I/O が律速になるのを防ぐ。
TMP_DIR = "/dev/shm/sand64_enum"


def load_prev_trails():
    """これまでに見つかった特性をすべて読み込む (重複は除く)"""
    seen = set()
    trails = []
    for path in sorted(glob.glob(os.path.join(PARTS_DIR, "*.txt"))):
        _, ts = parse(path)
        for t in ts:
            # 途中で kill された場合に備え, 13 状態そろったものだけを採用する
            if len(t) == NUM_ROUNDS + 1 and tuple(t) not in seen:
                seen.add(tuple(t))
                trails.append(t)
    return trails


def parse(path):
    wvec, trails, cur = None, [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# weight vector"):
                # 行末に別の情報が続く形式もあるので最初の "]" までを取る
                body = line[len("# weight vector"):]
                wvec = eval(body[:body.index("]") + 1].strip())
            elif line.startswith("# characteristic"):
                if cur:
                    trails.append(cur)
                cur = []
            elif line and not line.startswith("#"):
                p = line.split()
                if len(p) == 5:
                    cur.append((int(p[1] + p[2], 16), int(p[3] + p[4], 16)))
    if cur:
        trails.append(cur)
    return wvec, trails


def trail_bits(trail):
    """key_vars と同じ順序のビット列: dY^0 (32bit), そのあと dX^0..dX^R"""
    bits = [(trail[0][1] >> (31 - i)) & 1 for i in range(32)]
    for r in range(NUM_ROUNDS + 1):
        bits += [(trail[r][0] >> (31 - i)) & 1 for i in range(32)]
    return bits


def blocking_clause(trail, kv):
    return [(-v if b else v) for v, b in zip(kv, trail_bits(trail))]


def round1_active_set(trail):
    """ラウンド 1 で活性な SSb の位置集合 = dX^0 の非零ニブルの集合

    SSb の同時 DDT には din != 0 かつ重み 0 のエントリが存在しない
    (エントリ 16 は din = 0 のみ) ので, この集合は ssb_w_{j}_1 が真になる
    位置の集合と厳密に一致する。よってサブスライスの絞り込みは exact である。
    """
    dx0 = trail[0][0]
    S = set()
    for j in range(8):
        if any((dx0 >> b) & 1 for b in (31 - j, 23 - j, 15 - j, 7 - j)):
            S.add(j)
    return frozenset(S)


def job_key(idx, S):
    return f"v{idx:04d}_s" + "-".join(str(j) for j in S)


def enumerate_job(job):
    idx, vec, S, prev = job
    tag = job_key(idx, S)
    out_path = os.path.join(OUT_DIR, tag + ".txt")
    done_path = out_path + ".done"
    if os.path.exists(done_path):
        with open(done_path) as f:
            return tag, int(f.read().strip())

    t0 = time.time()
    model, variables = W.build()
    model.fix_weight_vector(vec)
    # ラウンド 1 の活性 SSb 位置を S に固定する
    for j in range(8):
        lit = variables[f"ssb_w_{j}_1"]
        model.clauses.append([lit] if j in S else [-lit])

    kv = key_vars(variables, NUM_ROUNDS)
    for t in prev:
        model.clauses.append(blocking_clause(t, kv))

    num_vars = max(max(abs(l) for l in c) for c in model.clauses)
    # ベース部分の DIMACS 文字列は一度だけ組み立て, 以後は追記分だけを連結する
    base_text = "".join(" ".join(map(str, c)) + " 0\n" for c in model.clauses)
    base_count = len(model.clauses)
    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"{tag}.cnf")

    found = 0
    block_text = []
    with open(out_path, "w") as out:
        out.write(f"# weight vector {vec}, round-1 active SSb positions {sorted(S)}\n")
        while True:
            with open(cnf, "w") as f:
                f.write(f"p cnf {num_vars} {base_count + len(block_text)}\n")
                f.write(base_text)
                f.write("".join(block_text))
            sat, pos = run_kissat(model.kissat_path, cnf)
            if not sat:
                break
            trail = extract_trail(pos, variables, NUM_ROUNDS)
            found += 1
            clause = [-v if v in pos else v for v in kv]
            block_text.append(" ".join(map(str, clause)) + " 0\n")
            out.write(f"\n# characteristic {found}\n")
            for r, (dx, dy) in enumerate(trail):
                out.write(f"{r:>3} {fmt(dx)} {fmt(dy)}\n")
            out.flush()
    if os.path.exists(cnf):
        os.remove(cnf)
    with open(done_path, "w") as f:
        f.write(str(found))
    print(f"{tag} -> {found} new ({time.time()-t0:.0f}s)", flush=True)
    return tag, found


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    with open(W.VEC_FILE) as f:
        vectors = json.load(f)
    prev = load_prev_trails()
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"resuming with {len(prev)} already-known characteristics as blocking clauses", flush=True)

    # サブスライスごとに, そこに属し得る既知特性だけをブロッキング節にする。
    # 他のスライスの特性は制約により既に除外されているので入れるだけ無駄であり,
    # CNF が数 MB 膨らんで探索が遅くなる原因になっていた。
    by_set = {}
    for t in prev:
        by_set.setdefault(round1_active_set(t), []).append(t)

    jobs = []
    for idx, vec in enumerate(vectors):
        size = vec[0] // 2            # ラウンド 1 の活性 SSb 個数 (重み 2 が size 個)
        for S in itertools.combinations(range(8), size):
            jobs.append((idx, vec, set(S), by_set.get(frozenset(S), [])))
    print(f"{len(jobs)} sub-slices, {workers} workers; "
          f"blocking clauses per slice: "
          f"{min(len(j[3]) for j in jobs)}..{max(len(j[3]) for j in jobs)}", flush=True)

    t0 = time.time()
    with mp.Pool(processes=workers) as pool:
        results = pool.map(enumerate_job, jobs, chunksize=1)
    new = sum(c for _, c in results)
    print("=" * 50)
    print(f"previously known : {len(prev)}")
    print(f"newly found      : {new}")
    print(f"TOTAL            : {len(prev) + new} characteristics of weight exactly 58 "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
