"""
SAND-128 の重み下限を確認するためのドライバ。

    python3 run_sand128_bound.py <rounds> <bound> <max_plaintext_ones> <max_output_ones>

<bound> は「総重み <= bound」の上限制約なので, UNSATISFIABLE であれば
その丸め条件のもとで重み <= bound の特性が存在しないことを意味する。
max_*_ones に 64 を与えると sequential_encoding_method() が m >= n で
早期 return するため, 平文側・暗号文側の絞り込みは無効になる（＝制約なし）。
"""
import importlib.util
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "SAND-128_DifferentialCryptanalysis_Kissat_Standalone_20260805.py")

spec = importlib.util.spec_from_file_location("sand128", SRC)
sand128 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sand128)

r = int(sys.argv[1])
b = int(sys.argv[2])
max_plaintext_ones = int(sys.argv[3])
max_output_ones = int(sys.argv[4])

model = sand128.SAND_128_SAT_Model()
print(f"SAND-128 forward search: rounds={r}, bound(total weight <=)={b}, "
      f"max_plaintext_ones={max_plaintext_ones}, max_output_ones={max_output_ones}")
model.verify_ssb_model()
print("-" * 50)

t0 = time.time()
model.run_without_key(r, b, max_plaintext_ones, max_output_ones)
print(f"elapsed: {time.time() - t0:.1f} s")
