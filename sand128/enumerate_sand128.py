"""
SAND-128 の r ラウンド・総重み <= b の差分特性を all-SAT で数え上げる。

    python3 enumerate_sand128.py <rounds> <bound> [出力ファイル]

SAND-128_DifferentialCryptanalysis_Kissat_Standalone_20260805.py の
順方向モデル (Search_Differential_Distingusher) をそのまま使い, 平文側・
暗号文側の 1 の個数を絞る制約は課さない。見つかった特性ごとに, その特性を
一意に定める変数への割り当てを否定する阻止節 (blocking clause) を追加して
再度解く, という操作を UNSATISFIABLE になるまで繰り返す。

特性を一意に定める変数:
    ΔX^0, ΔY^0 (x_i_0, 0 <= i < 128) と ΔX^r (x_i_r, 0 <= i < 64, 1 <= r <= R)
順方向モデルでは ΔY^r = ΔX^{r-1} なので, これで特性全体が決まる。

出力ファイルが既に存在する場合はそこから読み直して阻止節を張り直すので,
中断しても再開できる。
"""
import importlib.util
import os
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SAND-128_DifferentialCryptanalysis_Kissat_Standalone_20260805.py")

spec = importlib.util.spec_from_file_location("sand128", SRC)
sand128 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sand128)


def key_vars(variables, rounds):
    """特性を一意に定める変数の並び"""
    names = [f"x_{i}_0" for i in range(128)]
    for r in range(1, rounds + 1):
        names += [f"x_{i}_{r}" for i in range(64)]
    return [variables[n] for n in names]


def weight_vars(variables, rounds):
    return [variables[f"ssb_w_{i}_{r}"] for r in range(1, rounds + 1) for i in range(48)]


def write_cnf(path, clauses, num_vars):
    with open(path, "w") as f:
        f.write(f"p cnf {num_vars} {len(clauses)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in clauses))


def run_kissat(kissat, path):
    res = subprocess.run([kissat, "-q", path], capture_output=True, text=True)
    out = res.stdout
    if "s UNSATISFIABLE" in out:
        return False, None
    if "s SATISFIABLE" in out:
        positive = set()
        for line in out.split("\n"):
            if line.startswith("v "):
                for tok in line[2:].split():
                    v = int(tok)
                    if v > 0:
                        positive.add(v)
        return True, positive
    raise RuntimeError("unexpected kissat output:\n" + out + res.stderr)


def state_hex(model, variables, r, positive):
    dx = [1 if variables[f"x_{i}_{r}"] in positive else 0 for i in range(0, 64)]
    dy = [1 if variables[f"x_{i}_{r}"] in positive else 0 for i in range(64, 128)]
    return model.bits_to_hex(dx).replace(" ", ""), model.bits_to_hex(dy).replace(" ", "")


def parse_existing(path, rounds):
    """既存の出力ファイルから (特性ごとの ΔX/ΔY の16進表現) を読み戻す"""
    if not os.path.exists(path):
        return []
    trails, cur = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith("#"):
            if cur:
                trails.append(cur)
            cur = []
        elif line:
            parts = line.split()
            cur.append((int(parts[0]), parts[1], parts[2]))
    if cur:
        trails.append(cur)
    for t in trails:
        if len(t) != rounds + 1:
            raise RuntimeError(f"broken trail in {path}: {len(t)} lines, expected {rounds + 1}")
    return trails


def blocking_from_hex(trail, variables, rounds):
    """16進で読み戻した特性に対応する阻止節を作る"""
    clause = []
    for r, dx_hex, dy_hex in trail:
        value = int(dx_hex, 16)
        for i in range(64):
            v = variables[f"x_{i}_{r}"]
            clause.append(-v if (value >> (63 - i)) & 1 else v)
        if r == 0:
            value = int(dy_hex, 16)
            for i in range(64):
                v = variables[f"x_{i + 64}_0"]
                clause.append(-v if (value >> (63 - i)) & 1 else v)
    return clause


def main():
    rounds = int(sys.argv[1])
    bound = int(sys.argv[2])
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, f"trails_sand128_r{rounds}_w{bound}.txt")

    model = sand128.SAND_128_SAT_Model()
    model.verify_ssb_model()

    variables = model.generate_variables(rounds)
    model.Search_Differential_Distingusher(rounds, bound, variables)
    # 入力差分 != 0
    model.clauses.append([variables[f"x_{i}_0"] for i in range(128)])

    kvars = key_vars(variables, rounds)
    wvars = weight_vars(variables, rounds)
    num_vars = model.aux_var_counter
    base = list(model.clauses)

    print(f"SAND-128 enumeration: rounds={rounds}, total weight <= {bound}, no plaintext/ciphertext restriction")
    print(f"CNF: {num_vars} variables, {len(base)} clauses, {len(kvars)} key variables")
    print(f"output: {out_path}", flush=True)

    found = parse_existing(out_path, rounds)
    if found:
        print(f"resuming: {len(found)} characteristic(s) already in {out_path}", flush=True)
    blocking = [blocking_from_hex(t, variables, rounds) for t in found]

    weight_hist = {}
    cnf_path = tempfile.NamedTemporaryFile(suffix=".cnf", delete=False).name
    t_start = time.time()
    try:
        while True:
            write_cnf(cnf_path, base + blocking, num_vars)
            t0 = time.time()
            sat, positive = run_kissat(model.kissat_path, cnf_path)
            dt = time.time() - t0
            if not sat:
                total = time.time() - t_start
                print(f"UNSATISFIABLE after {len(found)} characteristic(s) ({dt:.1f} s, total {total:.1f} s)")
                print(f"=> exactly {len(found)} characteristic(s) with {rounds} rounds and total weight <= {bound}")
                for w in sorted(weight_hist):
                    print(f"   weight {w}: {weight_hist[w]}")
                return

            weight = sum(1 for v in wvars if v in positive)
            weight_hist[weight] = weight_hist.get(weight, 0) + 1
            trail = [(r,) + state_hex(model, variables, r, positive) for r in range(rounds + 1)]
            found.append(trail)

            with open(out_path, "a") as f:
                f.write(f"# characteristic {len(found)}  weight {weight}\n")
                for r, dx, dy in trail:
                    f.write(f"{r:>3} {dx} {dy}\n")

            blocking.append([(-v if v in positive else v) for v in kvars])
            print(f"[{len(found)}] weight {weight}  dX^0={trail[0][1]} dY^0={trail[0][2]} "
                  f"dX^{rounds}={trail[rounds][1]}  ({dt:.1f} s)", flush=True)
    finally:
        if os.path.exists(cnf_path):
            os.remove(cnf_path)


if __name__ == "__main__":
    main()
