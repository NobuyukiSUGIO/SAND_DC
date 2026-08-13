#!/usr/bin/env python3
"""
SAND-128 : Matsui 型の下界枝刈りを入れて r 段の最小重み B[r] を r = 1, 2, ... と順に求める。

枝刈りの原理
    R 段・総重み <= W の探索において, 先頭 k ラウンドの重みは
        （ラウンド 1..k の重み） <= W - B[R-k]
    を満たさなければならない。残る R-k ラウンドの重みは下界 B[R-k] 以上だからである。

実装がほぼ無料である理由
    元スクリプトの s_box_model() は self.vars にラウンド順で重みビットを追加する
    (ラウンド r のビットは添字 48(r-1)..48r-1)。したがって総重み <= W のために
    張られた sequential counter の補助変数 s[(48k, j)] が, そのまま
        「ラウンド 1..k の重みが j 以上」
    を意味する。よって上の制約は単位節
        [ -s[(48k, W - B[R-k] + 1)] ]
    を R-1 本追加するだけで表現でき, 新しい変数は一切要らない。

使い方
    python3 matsui_sand128.py <max_rounds> [timeout秒] [--no-matsui]

    --no-matsui : 枝刈りを外して同じ値を計算する (正しさの照合用)

B[r] は bounds_sand128.json に逐次保存され, 再実行時に読み直して続きから計算する。
"""
import importlib.util
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SAND-128_DifferentialCryptanalysis_Kissat_Standalone_20260805.py")
PER_ROUND = 48          # 1 ラウンドあたりの重みビット数 (16 SSb x 3)

spec = importlib.util.spec_from_file_location("sand128", SRC)
sand128 = importlib.util.module_from_spec(spec)
sys.modules["sand128"] = sand128
spec.loader.exec_module(sand128)


class CounterModel(sand128.SAND_128_SAT_Model):
    """sequential counter の補助変数を保持できるようにした派生クラス"""

    def __init__(self):
        super().__init__()
        self.counter_s = None
        self.counter_n = None

    def sequential_encoding_method(self, m):
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


def build(rounds, weight, bounds, use_matsui):
    """rounds 段・総重み <= weight のモデルを組み, Matsui 制約を加える"""
    m = CounterModel()
    variables = m.generate_variables(rounds)
    m.Search_Differential_Distingusher(rounds, weight, variables)
    # 状態差分は全 0 ではない (ΔX^0 = 0, ΔY^0 != 0 は許す)
    m.clauses.append([variables[f"x_{i}_0"] for i in range(128)])

    added = []
    if use_matsui:
        for k in range(1, rounds):
            b = bounds.get(rounds - k)
            if b is None:
                continue
            limit = weight - b                      # ラウンド 1..k の重みの上限
            j = limit + 1
            if limit < 0:
                # 残りだけで weight を超える -> この weight では解なし
                m.clauses.append([])                # 空節 = 即 UNSAT
                added.append((k, limit))
                break
            if (PER_ROUND * k, j) in m.counter_s:
                m.clauses.append([-m.counter_s[(PER_ROUND * k, j)]])
                added.append((k, limit))
    return m, variables, added


def run_kissat(kissat, path, timeout):
    try:
        res = subprocess.run([kissat, "-q", path], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None
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


TMP_DIR = "/dev/shm/sand128_matsui"


def test(rounds, weight, bounds, use_matsui, timeout):
    """rounds 段で総重み <= weight の特性が存在するか. (sat, 秒, 実際の重み, ラウンド重み)"""
    m, variables, added = build(rounds, weight, bounds, use_matsui)
    # 空節が入っていれば解かずに UNSAT
    if any(len(c) == 0 for c in m.clauses):
        return False, 0.0, None, None

    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"r{rounds}_w{weight}{'_nm' if not use_matsui else ''}.cnf")
    num_vars = m.aux_var_counter
    with open(cnf, "w") as f:
        f.write(f"p cnf {num_vars} {len(m.clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in m.clauses))

    t0 = time.time()
    sat, pos = run_kissat(m.kissat_path, cnf, timeout)
    dt = time.time() - t0
    if os.path.exists(cnf):
        os.remove(cnf)

    if sat is None:
        return None, dt, None, None
    if not sat:
        return False, dt, None, None

    rw = [sum(1 for i in range(PER_ROUND) if variables[f"ssb_w_{i}_{r}"] in pos)
          for r in range(1, rounds + 1)]
    return True, dt, sum(rw), rw


BOUND_FILE = os.environ.get("SAND128_BOUNDS", os.path.join(HERE, "bounds_sand128.json"))


def load_bounds():
    if os.path.exists(BOUND_FILE):
        return {int(k): v for k, v in json.load(open(BOUND_FILE)).items()}
    return {}


def save_bounds(bounds):
    with open(BOUND_FILE, "w") as f:
        json.dump({str(k): v for k, v in sorted(bounds.items())}, f, indent=1)


def find_bound(rounds, bounds, use_matsui, timeout):
    """B[rounds] を求める. B は単調非減少なので B[rounds-1] から上に探す"""
    lo = bounds.get(rounds - 1, 0)          # B[rounds] >= lo (確定した下界)
    calls = 0

    # まず SAT になる上界を指数的に探す (SAT 判定は概して安い)
    hi = None
    step = 0
    while True:
        w = lo + step
        sat, dt, tot, rw = test(rounds, w, bounds, use_matsui, timeout)
        calls += 1
        if sat is None:
            print(f"    r={rounds} w={w}: TIMEOUT ({dt:.0f}s) -> 中断", flush=True)
            return None, calls
        if sat:
            print(f"    r={rounds} w<={w}: SAT   ({dt:.1f}s) total={tot} rw={rw}", flush=True)
            hi = w
            break
        print(f"    r={rounds} w<={w}: UNSAT ({dt:.1f}s)", flush=True)
        lo = w + 1
        step = 1 if step == 0 else step * 2

    # (lo-1, hi] を二分する。lo-1 では UNSAT が確定している
    while lo < hi:
        mid = (lo + hi) // 2
        sat, dt, tot, rw = test(rounds, mid, bounds, use_matsui, timeout)
        calls += 1
        if sat is None:
            print(f"    r={rounds} w={mid}: TIMEOUT ({dt:.0f}s) -> B[{rounds}] in [{lo},{hi}]", flush=True)
            return None, calls
        if sat:
            print(f"    r={rounds} w<={mid}: SAT   ({dt:.1f}s) total={tot} rw={rw}", flush=True)
            hi = mid
        else:
            print(f"    r={rounds} w<={mid}: UNSAT ({dt:.1f}s)", flush=True)
            lo = mid + 1
    return hi, calls


def main():
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    timeout = float(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("-") else 3600.0
    use_matsui = "--no-matsui" not in sys.argv

    bounds = load_bounds()
    print(f"SAND-128 minimum-weight bounds, Matsui pruning: {'ON' if use_matsui else 'OFF'}, "
          f"per-call timeout {timeout:.0f}s")
    if bounds:
        print(f"loaded: {dict(sorted(bounds.items()))}")

    t_all = time.time()
    for r in range(1, max_rounds + 1):
        if r in bounds:
            continue
        print(f"  --- B[{r}] ---", flush=True)
        t0 = time.time()
        b, calls = find_bound(r, bounds, use_matsui, timeout)
        if b is None:
            print(f"  B[{r}] 未確定 (timeout) -> 打ち切り", flush=True)
            break
        bounds[r] = b
        save_bounds(bounds)
        print(f"  B[{r}] = {b}   ({calls} solver calls, {time.time()-t0:.1f}s)", flush=True)

    print("=" * 60)
    print(" r  B[r]  B[r]-B[r-1]")
    prev = None
    for r in sorted(bounds):
        d = "" if prev is None else str(bounds[r] - prev)
        print(f"{r:>3} {bounds[r]:>5}  {d:>5}")
        prev = bounds[r]
    print(f"total {time.time()-t_all:.1f}s; saved to {BOUND_FILE}")


if __name__ == "__main__":
    main()
