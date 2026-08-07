#!/usr/bin/env python3
"""
SAND-64 : 12段差分特性のうち重み (=-log2 p) がちょうど 58 のものを全列挙する。

使い方
    python3 enumerate_w58.py check <weight>   # 「重み <= weight の12段特性が存在するか」を判定
    python3 enumerate_w58.py enum  <weight>   # 「重み == weight」の12段特性を全列挙する

列挙の単位は「状態差分の系列」 (ΔX^0, ΔY^0), (ΔX^1, ΔY^1), ..., (ΔX^12, ΔY^12)
である。順方向モデルでは ΔY^r = ΔX^{r-1} が成り立つので, 系列は
    ΔY^0 と ΔX^0, ΔX^1, ..., ΔX^12   (32 * 14 = 448 bit)
で一意に決まる。この 448 変数に対してブロッキング節を追加しながら
kissat を繰り返し呼び出す (all-SAT)。

重み <= 58 の制約は元スクリプトの sequential counter をそのまま用い,
重み >= 58 の制約は同じ counter の補助変数 s[n][58]
(「先頭 n 個の重みビットのうち 58 個以上が真」) に単位節を課すことで実現する。
これにより「ちょうど 58」を追加変数なしで表現できる。
"""
import importlib.util
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


sand = _load("sand64_model", "SAND-64_20260729.py")
verifier = _load("sand64_verify", "verify_sand64_key_recovery.py")


class EnumModel(sand.SAND_64_SAT_Model):
    """sequential counter の補助変数を保持できるようにした派生クラス"""

    def __init__(self):
        super().__init__()
        self.counter_s = None
        self.counter_n = None

    def sequential_encoding_method(self, m):
        # 親クラスと同一のエンコードだが, 補助変数 s を self に保存する
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
        """重みの総和 >= k を課す (at-most-k のカウンタを共有する)"""
        s, n = self.counter_s, self.counter_n
        assert s is not None and (n, k) in s, "counter variable s[n][k] is missing"
        self.clauses.append([s[(n, k)]])


def build(num_rounds, weight, exact):
    m = EnumModel()
    variables = m.generate_variables(num_rounds)
    m.Search_Differential_Distingusher(num_rounds, weight, variables)
    # 入力差分は全 0 ではない
    m.clauses.append([variables[f"x_{i}_0"] for i in range(64)])
    if exact:
        m.add_at_least(weight)
    return m, variables


def key_vars(variables, num_rounds):
    """特性を一意に定める変数集合: ΔY^0 と ΔX^0..ΔX^R"""
    kv = [variables[f"x_{i}_0"] for i in range(32, 64)]
    for r in range(0, num_rounds + 1):
        kv += [variables[f"x_{i}_{r}"] for i in range(32)]
    return kv


def write_cnf(path, clauses, num_vars):
    with open(path, "w") as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in clauses))


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
    """(ΔX^r, ΔY^r) を 32bit 整数の組の列として取り出す"""
    trail = []
    for r in range(num_rounds + 1):
        dx = bits_to_word([1 if variables[f"x_{i}_{r}"] in pos else 0 for i in range(32)])
        dy = bits_to_word([1 if variables[f"x_{i}_{r}"] in pos else 0 for i in range(32, 64)])
        trail.append((dx, dy))
    return trail


def round_weights(pos, variables, num_rounds):
    return [sum(1 for i in range(24) if variables[f"ssb_w_{i}_{r}"] in pos)
            for r in range(1, num_rounds + 1)]


def fmt(word):
    return f"{word >> 16:04X} {word & 0xFFFF:04X}"


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "enum"
    weight = int(sys.argv[2]) if len(sys.argv) > 2 else 58
    num_rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 12

    exact = (mode == "enum")
    model, variables = build(num_rounds, weight, exact)
    num_vars = max(max(abs(l) for l in c) for c in model.clauses)
    base = model.clauses
    print(f"mode={mode}  rounds={num_rounds}  weight{'==' if exact else '<='}{weight}")
    print(f"variables: {num_vars}, clauses: {len(base)}", flush=True)

    cnf = os.path.join(HERE, f"_enum_{mode}_{weight}.cnf")
    kissat = model.kissat_path

    if mode == "check":
        write_cnf(cnf, base, num_vars)
        t0 = time.time()
        sat, pos = run_kissat(kissat, cnf)
        print(f"result: {'SATISFIABLE' if sat else 'UNSATISFIABLE'}  ({time.time()-t0:.1f} s)")
        if sat:
            tr = extract_trail(pos, variables, num_rounds)
            w = round_weights(pos, variables, num_rounds)
            for r, (dx, dy) in enumerate(tr):
                print(f"{r:>3} {fmt(dx)} {fmt(dy)}" + (f"  w={w[r-1]}" if r else ""))
            print("total weight =", sum(w))
        os.remove(cnf)
        return

    # ---- all-SAT enumeration ----
    kv = key_vars(variables, num_rounds)
    out_path = os.path.join(HERE, f"trails_w{weight}_r{num_rounds}.txt")
    found = []
    blocking = []
    t_start = time.time()
    with open(out_path, "w") as out:
        out.write(f"# SAND-64 {num_rounds}-round differential characteristics with weight exactly {weight}\n")
        while True:
            write_cnf(cnf, base + blocking, num_vars)
            t0 = time.time()
            sat, pos = run_kissat(kissat, cnf)
            dt = time.time() - t0
            if not sat:
                print(f"[{len(found)}] UNSAT after {dt:.1f}s -> enumeration complete", flush=True)
                break
            trail = extract_trail(pos, variables, num_rounds)
            w = round_weights(pos, variables, num_rounds)
            found.append((trail, w))
            blocking.append([-v if v in pos else v for v in kv])
            out.write(f"\n# characteristic {len(found)}  weights={w} sum={sum(w)}\n")
            for r, (dx, dy) in enumerate(trail):
                out.write(f"{r:>3} {fmt(dx)} {fmt(dy)}\n")
            out.flush()
            print(f"[{len(found)}] found in {dt:.1f}s  dX0={fmt(trail[0][0])} dY0={fmt(trail[0][1])} "
                  f"dX{num_rounds}={fmt(trail[-1][0])} dY{num_rounds}={fmt(trail[-1][1])} w={w}", flush=True)
    if os.path.exists(cnf):
        os.remove(cnf)
    print(f"total {len(found)} characteristics, {time.time()-t_start:.1f}s")
    print(f"written to {out_path}")


if __name__ == "__main__":
    main()
