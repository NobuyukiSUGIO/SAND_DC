#!/usr/bin/env python3
"""
SAND-64 : 重み (=-log2 p) がちょうど 58 の 12 段差分特性を全列挙する。

探索空間を「ラウンドごとの重みベクトル」 w = (w_1, ..., w_12), sum(w) = 58
で分割する。異なる w は互いに素で全体を覆うので, 各 w について独立に
all-SAT 列挙を行い和集合を取ればよい。w による分割は前半重みなどより
はるかに細かく, 16 コアの並列化が有効に効く。

    phase A : 実現可能な重みベクトル w をすべて列挙する (block on w)
    phase B : 各 w について特性を all-SAT 列挙する (block on 状態差分 448 bit)

使い方
    python3 enum_by_wvec.py vectors            # phase A
    python3 enum_by_wvec.py run [workers]      # phase B (phase A の結果を使う)
    python3 enum_by_wvec.py one <index>        # phase B を 1 つだけ実行 (デバッグ用)
"""
import json
import multiprocessing as mp
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from enumerate_w58 import (EnumModel, key_vars, write_cnf, run_kissat,
                           extract_trail, round_weights, fmt)

TOTAL_WEIGHT = 58
NUM_ROUNDS = 12
PER_ROUND_MAX = 24          # 1 ラウンドあたりの重みビット数 (8 SSb x 3)
VEC_FILE = os.path.join(HERE, "weight_vectors_58.json")
OUT_DIR = os.path.join(HERE, "trails_w58")


class WVecModel(EnumModel):
    """各ラウンドの重みに sequential counter を張り, 重みベクトルを直接扱えるようにする"""

    def full_counter(self, vars_list):
        """s[(i,j)] <-> 「先頭 i 個のうち j 個以上が真」を n = len(vars_list) 全域で構成する"""
        n = len(vars_list)
        s = {}
        for i in range(1, n + 1):
            for j in range(1, i + 1):
                self.aux_var_counter += 1
                s[(i, j)] = self.aux_var_counter

        self.clauses.append([-vars_list[0], s[(1, 1)]])
        self.clauses.append([vars_list[0], -s[(1, 1)]])

        for i in range(2, n + 1):
            # j = 1 : s[i][1] = s[i-1][1] OR x_i
            self.clauses.append([-s[(i - 1, 1)], s[(i, 1)]])
            self.clauses.append([-vars_list[i - 1], s[(i, 1)]])
            self.clauses.append([s[(i - 1, 1)], vars_list[i - 1], -s[(i, 1)]])
            for j in range(2, i + 1):
                if (i - 1, j) in s:
                    # s[i][j] = s[i-1][j] OR (s[i-1][j-1] AND x_i)
                    self.clauses.append([-s[(i - 1, j)], s[(i, j)]])
                    self.clauses.append([-s[(i - 1, j - 1)], -vars_list[i - 1], s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], s[(i - 1, j - 1)], -s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], vars_list[i - 1], -s[(i, j)]])
                else:
                    # j = i : s[i][i] = s[i-1][i-1] AND x_i
                    self.clauses.append([-s[(i - 1, j - 1)], -vars_list[i - 1], s[(i, j)]])
                    self.clauses.append([s[(i - 1, j - 1)], -s[(i, j)]])
                    self.clauses.append([vars_list[i - 1], -s[(i, j)]])
        return s

    def build_round_counters(self, variables):
        """ラウンド r の重み用カウンタ. self.round_s[r][v] <-> 「ラウンド r の重み >= v」"""
        self.round_s = {}
        for r in range(1, NUM_ROUNDS + 1):
            wv = [variables[f"ssb_w_{i}_{r}"] for i in range(PER_ROUND_MAX)]
            s = self.full_counter(wv)
            self.round_s[r] = {v: s[(PER_ROUND_MAX, v)] for v in range(1, PER_ROUND_MAX + 1)}

    def block_weight_vector(self, vec):
        """重みベクトル vec をひとつ禁止する節を追加する"""
        clause = []
        for r, v in enumerate(vec, start=1):
            if v >= 1:
                clause.append(-self.round_s[r][v])
            if v + 1 <= PER_ROUND_MAX:
                clause.append(self.round_s[r][v + 1])
        self.clauses.append(clause)

    def fix_weight_vector(self, vec):
        """重みベクトルを vec に固定する単位節を追加する"""
        for r, v in enumerate(vec, start=1):
            if v >= 1:
                self.clauses.append([self.round_s[r][v]])
            if v + 1 <= PER_ROUND_MAX:
                self.clauses.append([-self.round_s[r][v + 1]])


def build(with_counters=True):
    m = WVecModel()
    variables = m.generate_variables(NUM_ROUNDS)
    m.Search_Differential_Distingusher(NUM_ROUNDS, TOTAL_WEIGHT, variables)
    m.clauses.append([variables[f"x_{i}_0"] for i in range(64)])
    m.add_at_least(TOTAL_WEIGHT)
    if with_counters:
        m.build_round_counters(variables)
    return m, variables


# --------------------------------------------------------------------------
# phase A : 実現可能な重みベクトルの列挙
# --------------------------------------------------------------------------
def phase_a_slice(a):
    """ラウンド 1 の重みが a である重みベクトルをすべて列挙する"""
    t0 = time.time()
    model, variables = build()
    # w_1 == a に固定
    if a >= 1:
        model.clauses.append([model.round_s[1][a]])
    if a + 1 <= PER_ROUND_MAX:
        model.clauses.append([-model.round_s[1][a + 1]])
    num_vars = max(max(abs(l) for l in c) for c in model.clauses)
    cnf = os.path.join(HERE, f"_wvec_{a}.cnf")
    vectors = []
    while True:
        write_cnf(cnf, model.clauses, num_vars)
        sat, pos = run_kissat(model.kissat_path, cnf)
        if not sat:
            break
        vec = round_weights(pos, variables, NUM_ROUNDS)
        assert sum(vec) == TOTAL_WEIGHT and vec[0] == a, vec
        vectors.append(vec)
        model.block_weight_vector(vec)
    if os.path.exists(cnf):
        os.remove(cnf)
    print(f"  w1={a:>2}: {len(vectors)} vectors ({time.time()-t0:.0f}s)", flush=True)
    return vectors


def phase_a(workers=15):
    t0 = time.time()
    print(f"phase A: enumerating realizable weight vectors, {workers} workers", flush=True)
    with mp.Pool(processes=workers) as pool:
        chunks = pool.map(phase_a_slice, list(range(0, PER_ROUND_MAX + 1)), chunksize=1)
    vectors = [v for c in chunks for v in c]
    with open(VEC_FILE, "w") as f:
        json.dump(vectors, f)
    print(f"phase A done: {len(vectors)} realizable weight vectors, {time.time()-t0:.0f}s")
    return vectors


# --------------------------------------------------------------------------
# phase B : 重みベクトルごとの特性列挙
# --------------------------------------------------------------------------
def enumerate_vector(job):
    idx, vec = job
    out_path = os.path.join(OUT_DIR, f"v{idx:04d}.txt")
    done_path = out_path + ".done"
    if os.path.exists(done_path):
        with open(done_path) as f:
            return idx, int(f.read().strip())

    t0 = time.time()
    model, variables = build()
    model.fix_weight_vector(vec)
    num_vars = max(max(abs(l) for l in c) for c in model.clauses)
    base = model.clauses
    kv = key_vars(variables, NUM_ROUNDS)
    cnf = os.path.join(HERE, f"_v{idx:04d}.cnf")

    found = 0
    blocking = []
    with open(out_path, "w") as out:
        out.write(f"# weight vector {vec} (sum={sum(vec)})\n")
        while True:
            write_cnf(cnf, base + blocking, num_vars)
            sat, pos = run_kissat(model.kissat_path, cnf)
            if not sat:
                break
            trail = extract_trail(pos, variables, NUM_ROUNDS)
            found += 1
            blocking.append([-v if v in pos else v for v in kv])
            out.write(f"\n# characteristic {found}\n")
            for r, (dx, dy) in enumerate(trail):
                out.write(f"{r:>3} {fmt(dx)} {fmt(dy)}\n")
            out.flush()
    if os.path.exists(cnf):
        os.remove(cnf)
    with open(done_path, "w") as f:
        f.write(str(found))
    print(f"v{idx:04d} {vec} -> {found} ({time.time()-t0:.0f}s)", flush=True)
    return idx, found


def phase_b(workers):
    with open(VEC_FILE) as f:
        vectors = json.load(f)
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = list(enumerate(vectors))
    # 重い仕事を先に投げるためのヒューリスティック: 重みの分散が小さいものほど自由度が高い
    jobs.sort(key=lambda j: -max(j[1]))
    print(f"phase B: {len(jobs)} weight vectors, {workers} workers", flush=True)
    t0 = time.time()
    with mp.Pool(processes=workers) as pool:
        results = pool.map(enumerate_vector, jobs, chunksize=1)
    total = sum(c for _, c in results)
    print("=" * 50)
    print(f"TOTAL = {total} characteristics of weight exactly {TOTAL_WEIGHT} "
          f"({time.time()-t0:.0f}s)")
    return total


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "vectors"
    if cmd == "vectors":
        phase_a()
    elif cmd == "run":
        phase_b(int(sys.argv[2]) if len(sys.argv) > 2 else 15)
    elif cmd == "one":
        with open(VEC_FILE) as f:
            vectors = json.load(f)
        i = int(sys.argv[2])
        os.makedirs(OUT_DIR, exist_ok=True)
        enumerate_vector((i, vectors[i]))


if __name__ == "__main__":
    main()
