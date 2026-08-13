#!/usr/bin/env python3
"""CNF フィンガープリントによる回帰テスト。

各モデル構築関数が生成する CNF を SHA-256 で指紋化し, tests/cnf_baseline.json
に記録した値と比較する。CNF の生成は完全に決定的なので, 指紋が一致すれば
探索の挙動が変わっていないことを SAT を解かずに確認できる。

リファクタリング (例: sand64/ と sand128/ の共通部分を common/ に切り出す)
の前後でこれを走らせ, 差が出なければ振る舞いは保存されている。

    python3 tests/cnf_fingerprint.py            # 検査 (数秒)
    python3 tests/cnf_fingerprint.py --update   # 基準値を作り直す
    python3 tests/cnf_fingerprint.py --verbose  # 各ケースの内訳も表示

基準値を作り直すのは, CNF の変化が意図したものだと確認できたときだけにする。
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
BASELINE = os.path.join(HERE, "cnf_baseline.json")

# 各ケース: (識別子, 変種ディレクトリ, モジュール, 呼び出し式)
# 呼び出し式のなかで M はそのモジュール, BOUNDS はその変種の bounds を指す。
CASES = [
    # --- Matsui 枝刈り付きの下界探索 (matsui_sand*.py) ---
    ("sand64/matsui/r4_w8",        "sand64",  "matsui_sand64",
     "M.build(4, 8, BOUNDS, True)"),
    ("sand64/matsui/r6_w18",       "sand64",  "matsui_sand64",
     "M.build(6, 18, BOUNDS, True)"),
    ("sand64/matsui/r8_w32",       "sand64",  "matsui_sand64",
     "M.build(8, 32, BOUNDS, True)"),
    ("sand64/matsui/r8_w32_nom",   "sand64",  "matsui_sand64",
     "M.build(8, 32, BOUNDS, False)"),
    ("sand128/matsui/r4_w8",       "sand128", "matsui_sand128",
     "M.build(4, 8, BOUNDS, True)"),
    ("sand128/matsui/r6_w18",      "sand128", "matsui_sand128",
     "M.build(6, 18, BOUNDS, True)"),
    ("sand128/matsui/r8_w32",      "sand128", "matsui_sand128",
     "M.build(8, 32, BOUNDS, True)"),
    ("sand128/matsui/r8_w32_nom",  "sand128", "matsui_sand128",
     "M.build(8, 32, BOUNDS, False)"),

    # --- ちょうど k の基数制約を使う列挙モデル (enumerate_w58.py) ---
    ("sand64/enum/r6_w18_exact",   "sand64",  "enumerate_w58",
     "M.build(6, 18, True)"),
    ("sand64/enum/r6_w18_le",      "sand64",  "enumerate_w58",
     "M.build(6, 18, False)"),

    # --- 重みベクトルで分割する列挙モデル (enum_by_wvec*.py) ---
    ("sand64/wvec/r12_w58",        "sand64",  "enum_by_wvec",
     "M.build(True)"),
    ("sand128/wvec/r19_w120",      "sand128", "enum_by_wvec_sand128",
     "M.build(True)"),

    # --- 鍵回復用の境界制約 (keyrecovery_sand128.py) ---
    # 平文側と暗号文側で別の端を縛るので, 取り違えると指紋が変わる。
    ("sand128/kr/r6_w18_free",     "sand128", "keyrecovery_sand128",
     "M.build(6, 18, BOUNDS)"),
    ("sand128/kr/r6_w18_pt4",      "sand128", "keyrecovery_sand128",
     "M.build(6, 18, BOUNDS, max_nib_pt=4)"),
    ("sand128/kr/r6_w18_ct4",      "sand128", "keyrecovery_sand128",
     "M.build(6, 18, BOUNDS, max_nib_ct=4)"),
    ("sand128/kr/r6_w18_pt4_ct6",  "sand128", "keyrecovery_sand128",
     "M.build(6, 18, BOUNDS, max_nib_pt=4, max_nib_ct=6)"),
]

# bounds をどのモジュールから読むか (keyrecovery は matsui のものを使う)
BOUNDS_MODULE = {"sand64": "matsui_sand64", "sand128": "matsui_sand128"}

WORKER = r'''
import hashlib, importlib, json, os, sys
subdir, module, expr, bounds_module = sys.argv[1:5]
sys.path.insert(0, os.path.join(%r, subdir))
M = importlib.import_module(module)
BOUNDS = importlib.import_module(bounds_module).load_bounds()
res = eval(expr)
m = res[0] if isinstance(res, tuple) else res
h = hashlib.sha256()
h.update(f"{m.aux_var_counter}\n".encode())
for c in m.clauses:
    h.update((" ".join(map(str, c)) + "\n").encode())
print(json.dumps({"vars": m.aux_var_counter, "clauses": len(m.clauses),
                  "sha256": h.hexdigest()}))
''' % (REPO,)


def fingerprint(case):
    """1 ケースを別プロセスで構築する。

    変種ごとのモジュールが同じ論理名で標準モデルを sys.modules に登録するため,
    同一プロセスで両方を読み込むと互いを汚染する。プロセスを分けて回避する。
    """
    ident, subdir, module, expr = case
    out = subprocess.run(
        [sys.executable, "-c", WORKER, subdir, module, expr, BOUNDS_MODULE[subdir]],
        capture_output=True, text=True, cwd=REPO)
    if out.returncode != 0:
        return {"error": (out.stderr.strip().splitlines() or ["(no output)"])[-1]}
    return json.loads(out.stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="基準値を作り直す")
    ap.add_argument("--verbose", action="store_true", help="一致したケースも表示する")
    args = ap.parse_args()

    baseline = {}
    if not args.update:
        if not os.path.exists(BASELINE):
            print(f"基準値がありません: {BASELINE}\n"
                  f"  python3 tests/cnf_fingerprint.py --update  で作成してください")
            return 2
        with open(BASELINE) as f:
            baseline = json.load(f)["cases"]

    results = {}
    failed = []
    for case in CASES:
        ident = case[0]
        got = fingerprint(case)
        results[ident] = got
        if args.update:
            mark, note = "rec ", ""
        elif "error" in got:
            mark, note = "FAIL", got["error"]
            failed.append(ident)
        elif ident not in baseline:
            mark, note = "NEW ", "基準値になし (--update で登録)"
            failed.append(ident)
        elif got == baseline[ident]:
            mark, note = "ok  ", ""
        else:
            b = baseline[ident]
            mark = "FAIL"
            note = (f"vars {b['vars']}->{got['vars']}, "
                    f"clauses {b['clauses']}->{got['clauses']}, "
                    f"sha256 {b['sha256'][:12]}->{got['sha256'][:12]}")
            failed.append(ident)
        if args.verbose or mark != "ok  ":
            extra = f"  {note}" if note else ""
            counts = ("" if "error" in got
                      else f"  vars={got['vars']:>7} clauses={got['clauses']:>8}"
                           f" {got['sha256'][:16]}")
            print(f"  [{mark}] {ident:<28}{counts}{extra}")

    if args.update:
        with open(BASELINE, "w") as f:
            json.dump({"note": "python3 tests/cnf_fingerprint.py で検査する。"
                               "CNF 生成は決定的なので, 指紋が一致すれば探索の"
                               "挙動は変わっていない。",
                       "cases": results}, f, indent=1, sort_keys=True)
            f.write("\n")
        print(f"\n{len(results)} ケースの基準値を {BASELINE} に記録しました")
        return 0

    print()
    if failed:
        print(f"{len(failed)}/{len(CASES)} ケースが不一致: " + ", ".join(failed))
        return 1
    print(f"{len(CASES)}/{len(CASES)} ケース一致 — CNF は変化していません")
    return 0


if __name__ == "__main__":
    sys.exit(main())
