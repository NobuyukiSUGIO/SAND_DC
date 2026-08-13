#!/usr/bin/env python3
"""
SAND-128 : 重み (=-log2 p) がちょうど W の R 段差分特性を全列挙する。

SAND_DC (SAND-64) の enum_by_wvec.py と同じ戦略を SAND-128 に移植したもの。
探索空間を「ラウンドごとの重みベクトル」 w = (w_1, ..., w_R), sum(w) = W で
分割する。異なる w は互いに素で全体を覆うので, 各 w について独立に all-SAT
列挙を行い和集合を取ればよい。

    phase A : 実現可能な重みベクトル w をすべて列挙する (block on w)
    phase B : 各 w について特性を all-SAT 列挙する (block on 状態差分)

使い方
    python3 enum_by_wvec_sand128.py vectors [workers]   # phase A
    python3 enum_by_wvec_sand128.py run     [workers]   # phase B
    python3 enum_by_wvec_sand128.py one     <index>     # phase B を 1 つだけ (デバッグ用)

ラウンド数と総重みは環境変数で変更できる (既定は R=19, W=120):
    SAND128_ROUNDS=2 SAND128_WEIGHT=2 python3 enum_by_wvec_sand128.py vectors
"""
import importlib.util
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SAND-128_DifferentialCryptanalysis_Kissat_Standalone_20260805.py")

NUM_ROUNDS = int(os.environ.get("SAND128_ROUNDS", "19"))
TOTAL_WEIGHT = int(os.environ.get("SAND128_WEIGHT", "120"))
# SAND128_EXACT=1 : 総重みちょうど TOTAL_WEIGHT / =0 : 総重み <= TOTAL_WEIGHT
EXACT = os.environ.get("SAND128_EXACT", "1") != "0"
PER_ROUND_MAX = 48          # 1 ラウンドあたりの重みビット数 (16 SSb x 3)
NUM_SSB = 16                # 1 ラウンドあたりの SSb 個数

TAG = f"r{NUM_ROUNDS}_w{'' if EXACT else 'le'}{TOTAL_WEIGHT}"
VEC_FILE = os.path.join(HERE, f"weight_vectors_sand128_{TAG}.json")
OUT_DIR = os.path.join(HERE, f"trails_sand128_{TAG}")
# CNF は毎回書き直すので tmpfs (RAM) に置き, ディスク I/O が律速になるのを防ぐ
TMP_DIR = f"/dev/shm/sand128_enum_{TAG}"

spec = importlib.util.spec_from_file_location("sand128", SRC)
sand128 = importlib.util.module_from_spec(spec)
sys.modules["sand128"] = sand128
spec.loader.exec_module(sand128)


# --------------------------------------------------------------------------
# モデル : sequential counter の補助変数を保持して exactly-k / 重みベクトルを扱う
# --------------------------------------------------------------------------
class WVecModel(sand128.SAND_128_SAT_Model):
    def __init__(self):
        super().__init__()
        self.counter_s = None
        self.counter_n = None
        self.round_s = None

    def sequential_encoding_method(self, m):
        """親クラスと同一のエンコードだが, 補助変数 s を self に保存する"""
        n = len(self.vars)
        if m >= n:
            self.counter_s, self.counter_n = {}, n
            return

        s = {}
        for i in range(1, n + 1):
            for j in range(1, min(i, m + 1) + 1):
                self.aux_var_counter += 1
                s[(i, j)] = self.aux_var_counter

        self.clauses.append([-self.vars[0], s[(1, 1)]])
        self.clauses.append([self.vars[0], -s[(1, 1)]])
        for j in range(2, min(2, m + 1) + 1):
            if (1, j) in s:
                self.clauses.append([-s[(1, j)]])

        for i in range(2, n + 1):
            if (i - 1, 1) in s and (i, 1) in s:
                self.clauses.append([-s[(i - 1, 1)], s[(i, 1)]])
                self.clauses.append([-self.vars[i - 1], s[(i, 1)]])
                self.clauses.append([s[(i - 1, 1)], self.vars[i - 1], -s[(i, 1)]])
            for j in range(2, min(i, m + 1) + 1):
                if (i - 1, j) in s and (i, j) in s and (i - 1, j - 1) in s:
                    self.clauses.append([-s[(i - 1, j)], s[(i, j)]])
                    self.clauses.append([-s[(i - 1, j - 1)], -self.vars[i - 1], s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], s[(i - 1, j - 1)], -s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], self.vars[i - 1], -s[(i, j)]])
                elif (i, j) in s and (i - 1, j - 1) in s:
                    self.clauses.append([-s[(i - 1, j - 1)], -self.vars[i - 1], s[(i, j)]])
                    self.clauses.append([s[(i - 1, j - 1)], -s[(i, j)]])
                    self.clauses.append([self.vars[i - 1], -s[(i, j)]])

        if m + 1 <= n and (n, m + 1) in s:
            self.clauses.append([-s[(n, m + 1)]])

        self.counter_s, self.counter_n = s, n

    def add_at_least(self, k):
        """重みの総和 >= k を課す (at-most-k のカウンタを共有するので追加変数は不要)"""
        s, n = self.counter_s, self.counter_n
        assert s is not None and (n, k) in s, "counter variable s[n][k] is missing"
        self.clauses.append([s[(n, k)]])

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
            self.clauses.append([-s[(i - 1, 1)], s[(i, 1)]])
            self.clauses.append([-vars_list[i - 1], s[(i, 1)]])
            self.clauses.append([s[(i - 1, 1)], vars_list[i - 1], -s[(i, 1)]])
            for j in range(2, i + 1):
                if (i - 1, j) in s:
                    self.clauses.append([-s[(i - 1, j)], s[(i, j)]])
                    self.clauses.append([-s[(i - 1, j - 1)], -vars_list[i - 1], s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], s[(i - 1, j - 1)], -s[(i, j)]])
                    self.clauses.append([s[(i - 1, j)], vars_list[i - 1], -s[(i, j)]])
                else:
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
    # 入力差分は全 0 ではない
    m.clauses.append([variables[f"x_{i}_0"] for i in range(128)])
    if EXACT:
        m.add_at_least(TOTAL_WEIGHT)      # <= と合わせて「ちょうど TOTAL_WEIGHT」
    if with_counters:
        m.build_round_counters(variables)
    return m, variables


# --------------------------------------------------------------------------
# 共通ユーティリティ
# --------------------------------------------------------------------------
def key_vars(variables, num_rounds):
    """特性を一意に定める変数集合: ΔX^0, ΔY^0 と ΔX^1..ΔX^R

    順方向モデルでは ΔY^r = ΔX^{r-1} なので, これで特性全体が決まる。
    """
    kv = [variables[f"x_{i}_0"] for i in range(128)]
    for r in range(1, num_rounds + 1):
        kv += [variables[f"x_{i}_{r}"] for i in range(64)]
    return kv


def run_kissat(kissat, path):
    res = subprocess.run([kissat, "-q", path], capture_output=True, text=True)
    out = res.stdout
    if "s SATISFIABLE" in out:
        pos = set()
        for line in out.split("\n"):
            if line.startswith("v "):
                for tok in line[2:].split():
                    v = int(tok)
                    if v > 0:
                        pos.add(v)
        return True, pos
    if "s UNSATISFIABLE" in out:
        return False, None
    raise RuntimeError("unexpected kissat output:\n" + out + res.stderr)


def bits_to_word(bits):
    v = 0
    for b in bits:
        v = (v << 1) | b
    return v


def extract_trail(pos, variables, num_rounds):
    """(ΔX^r, ΔY^r) を 64bit 整数の組の列として取り出す"""
    trail = []
    for r in range(num_rounds + 1):
        dx = bits_to_word([1 if variables[f"x_{i}_{r}"] in pos else 0 for i in range(0, 64)])
        dy = bits_to_word([1 if variables[f"x_{i}_{r}"] in pos else 0 for i in range(64, 128)])
        trail.append((dx, dy))
    return trail


def round_weights(pos, variables, num_rounds):
    return [sum(1 for i in range(PER_ROUND_MAX) if variables[f"ssb_w_{i}_{r}"] in pos)
            for r in range(1, num_rounds + 1)]


def fmt(word):
    return f"{word:016X}"


def write_cnf(path, clauses, num_vars):
    with open(path, "w") as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in clauses))


# --------------------------------------------------------------------------
# phase A : 実現可能な重みベクトルの列挙
# --------------------------------------------------------------------------
def phase_a_slice(a):
    """ラウンド 1 の重みが a である重みベクトルをすべて列挙する"""
    t0 = time.time()
    model, variables = build()
    if a >= 1:
        model.clauses.append([model.round_s[1][a]])
    if a + 1 <= PER_ROUND_MAX:
        model.clauses.append([-model.round_s[1][a + 1]])
    num_vars = model.aux_var_counter
    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"_wvec_{a}.cnf")
    vectors = []
    while True:
        write_cnf(cnf, model.clauses, num_vars)
        sat, pos = run_kissat(model.kissat_path, cnf)
        if not sat:
            break
        vec = round_weights(pos, variables, NUM_ROUNDS)
        if EXACT:
            assert sum(vec) == TOTAL_WEIGHT and vec[0] == a, vec
        else:
            assert sum(vec) <= TOTAL_WEIGHT and vec[0] == a, vec
        vectors.append(vec)
        model.block_weight_vector(vec)
    if os.path.exists(cnf):
        os.remove(cnf)
    print(f"  w1={a:>2}: {len(vectors)} vectors ({time.time()-t0:.0f}s)", flush=True)
    return vectors


def phase_a(workers):
    t0 = time.time()
    print(f"phase A: R={NUM_ROUNDS}, total weight {'==' if EXACT else '<='} {TOTAL_WEIGHT}, "
          f"{workers} workers", flush=True)
    with mp.Pool(processes=workers) as pool:
        chunks = pool.map(phase_a_slice, list(range(0, PER_ROUND_MAX + 1)), chunksize=1)
    vectors = [v for c in chunks for v in c]
    with open(VEC_FILE, "w") as f:
        json.dump(vectors, f)
    print(f"phase A done: {len(vectors)} realizable weight vectors, {time.time()-t0:.0f}s")
    print(f"written to {VEC_FILE}")
    return vectors


# --------------------------------------------------------------------------
# phase B : 重みベクトルごとの特性列挙
# --------------------------------------------------------------------------
def enumerate_vector(job):
    idx, vec = job
    out_path = os.path.join(OUT_DIR, f"v{idx:05d}.txt")
    done_path = out_path + ".done"
    if os.path.exists(done_path):
        with open(done_path) as f:
            return idx, int(f.read().strip())

    t0 = time.time()
    model, variables = build()
    model.fix_weight_vector(vec)
    num_vars = model.aux_var_counter
    kv = key_vars(variables, NUM_ROUNDS)

    # ベース部分の DIMACS 文字列は一度だけ組み立て, 以後は追記分だけを連結する
    base_text = "".join(" ".join(map(str, c)) + " 0\n" for c in model.clauses)
    base_count = len(model.clauses)
    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"v{idx:05d}.cnf")

    found = 0
    block_text = []
    with open(out_path, "w") as out:
        out.write(f"# weight vector {vec} (sum={sum(vec)})\n")
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
    print(f"v{idx:05d} {vec} -> {found} ({time.time()-t0:.0f}s)", flush=True)
    return idx, found


def phase_b(workers):
    with open(VEC_FILE) as f:
        vectors = json.load(f)
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = list(enumerate(vectors))
    # 重い仕事を先に投げるヒューリスティック: 最大ラウンド重みが大きいものほど自由度が高い
    jobs.sort(key=lambda j: -max(j[1]))
    print(f"phase B: {len(jobs)} weight vectors, {workers} workers", flush=True)
    t0 = time.time()
    with mp.Pool(processes=workers) as pool:
        results = pool.map(enumerate_vector, jobs, chunksize=1)
    total = sum(c for _, c in results)
    print("=" * 50)
    print(f"TOTAL = {total} characteristics of weight {'exactly' if EXACT else 'at most'} "
          f"{TOTAL_WEIGHT} over {NUM_ROUNDS} rounds ({time.time()-t0:.0f}s)")
    return total


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "vectors"
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    if cmd == "vectors":
        phase_a(workers)
    elif cmd == "run":
        phase_b(workers)
    elif cmd == "one":
        with open(VEC_FILE) as f:
            vectors = json.load(f)
        i = int(sys.argv[2])
        os.makedirs(OUT_DIR, exist_ok=True)
        enumerate_vector((i, vectors[i]))
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
