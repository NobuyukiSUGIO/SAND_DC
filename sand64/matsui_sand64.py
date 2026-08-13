#!/usr/bin/env python3
"""
SAND-64 : Matsui 型の下界枝刈りを入れて r 段の最小重み B[r] を r = 1, 2, ... と順に求める。

matsui_sand128.py と同じ手法を SAND-64 に適用したもの。SAND-64 は 8 個の SSb を
持つので, 重みビットは 1 ラウンドあたり 8 x 3 = 24 個であり, 逐次カウンタの
補助変数 s[(24k, j)] が「ラウンド 1..k の重みが j 以上」を意味する。

使い方
    python3 matsui_sand64.py <max_rounds> [timeout秒] [--no-matsui]

B[r] は bounds_sand64.json に逐次保存され, 再実行時に読み直して続きから計算する。
"""
import importlib.util
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SAND-64_20260729.py")
PER_ROUND = 24          # 1 ラウンドあたりの重みビット数 (8 SSb x 3)

spec = importlib.util.spec_from_file_location("sand64", SRC)
sand64 = importlib.util.module_from_spec(spec)
sys.modules["sand64"] = sand64
spec.loader.exec_module(sand64)


class CounterModel(sand64.SAND_64_SAT_Model):
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
    m = CounterModel()
    variables = m.generate_variables(rounds)
    m.Search_Differential_Distingusher(rounds, weight, variables)
    # 状態差分は全 0 ではない
    m.clauses.append([variables[f"x_{i}_0"] for i in range(64)])

    added = []
    if use_matsui:
        for k in range(1, rounds):
            b = bounds.get(rounds - k)
            if b is None:
                continue
            limit = weight - b
            j = limit + 1
            if limit < 0:
                m.clauses.append([])          # 空節 = 即 UNSAT
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


TMP_DIR = "/dev/shm/sand64_matsui"


def test(rounds, weight, bounds, use_matsui, timeout):
    m, variables, added = build(rounds, weight, bounds, use_matsui)
    if any(len(c) == 0 for c in m.clauses):
        return False, 0.0, None, None

    os.makedirs(TMP_DIR, exist_ok=True)
    cnf = os.path.join(TMP_DIR, f"r{rounds}_w{weight}{'_nm' if not use_matsui else ''}.cnf")
    with open(cnf, "w") as f:
        f.write(f"p cnf {m.aux_var_counter} {len(m.clauses)}\n")
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


BOUND_FILE = os.environ.get("SAND64_BOUNDS", os.path.join(HERE, "bounds_sand64.json"))


def load_bounds():
    if os.path.exists(BOUND_FILE):
        return {int(k): v for k, v in json.load(open(BOUND_FILE)).items()}
    return {}


def save_bounds(bounds):
    with open(BOUND_FILE, "w") as f:
        json.dump({str(k): v for k, v in sorted(bounds.items())}, f, indent=1)


def find_bound(rounds, bounds, use_matsui, timeout):
    lo = bounds.get(rounds - 1, 0)
    calls = 0
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
    max_rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    timeout = float(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("-") else 3600.0
    use_matsui = "--no-matsui" not in sys.argv

    bounds = load_bounds()
    print(f"SAND-64 minimum-weight bounds, Matsui pruning: {'ON' if use_matsui else 'OFF'}, "
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
