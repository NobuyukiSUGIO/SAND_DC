#!/usr/bin/env python3
"""
選択暗号文攻撃で使える鍵非依存のフィルタ条件を, 切詰差分の後方伝搬で求める。

論文 Fig.4 の「赤 (活性) / 灰 (未知)」の伝搬を 3 値論理 (0 / 1 / 未知) で計算する。

要点: 入力差分 delta が非零でも, 導関数 D_delta SSb の出力ビットの一部は定数に
なる。例えば delta = 0x4 では
    D_delta N0(u) = (u3u1, 1, 0, u3),  D_delta N1(u) = (0, 1, u0, 0)
なので 8 ビット中 5 ビットが確定する。これが論文の 2^-8, 2^-23 というフィルタの
源であり, ニブル単位で「非零なら全部未知」と扱うと強さを大幅に見誤る。

    python3 filters_sand128.py validate          # 論文 (SAND-64) の数値で検証
    python3 filters_sand128.py sand128 <dX_s> <dY_s> [steps]
"""
import sys

import keyrecovery_sand128 as K

UNK = None


def derivative_pattern(delta):
    """
    入力差分 delta に対する D_delta SSb の 8 出力ビットの確定/未知パターン。
    返り値 pat[0..3] = G0 の y3..y0, pat[4..7] = G1 の y3..y0 (0 / 1 / None)
    """
    if delta == 0:
        return [0] * 8
    vals = [K.ssb_value(u) ^ K.ssb_value(u ^ delta) for u in range(16)]
    pat = []
    for b in range(7, -1, -1):
        s = {(v >> b) & 1 for v in vals}
        pat.append(s.pop() if len(s) == 1 else UNK)
    return pat


def word_to_arr(word, nbits):
    """64 (32) ビット整数を index 0 = MSB の配列にする"""
    return [(word >> (nbits - 1 - i)) & 1 for i in range(nbits)]


def arr_to_str(arr, ncols):
    """行ごとに 0/1/? を並べた文字列 (上から u3 行, u2 行, u1 行, u0 行)"""
    rows = []
    for k in range(4):
        rows.append("".join("?" if arr[k * ncols + j] is UNK else str(arr[k * ncols + j])
                            for j in range(ncols)))
    return " ".join(rows)


def known_mask(arr, ncols):
    """確定しているビットのマスクと, そのうち 1 のビットのマスク"""
    nbits = 4 * ncols
    m = 0
    v = 0
    for i in range(nbits):
        if arr[i] is not UNK:
            m |= 1 << (nbits - 1 - i)
            if arr[i] == 1:
                v |= 1 << (nbits - 1 - i)
    return m, v


def possible_deltas(arr, j, ncols):
    """ニブル j の既知ビットと整合する入力差分 delta の全体"""
    d = 0
    unk = []
    for k, b in enumerate((3, 2, 1, 0)):
        x = arr[k * ncols + j]
        if x is UNK:
            unk.append(b)
        elif x:
            d |= 1 << b
    out = []
    for code in range(1 << len(unk)):
        dd = d
        for t, b in enumerate(unk):
            if (code >> t) & 1:
                dd |= 1 << b
        out.append(dd)
    return out


def combined_pattern(deltas):
    """
    取り得る delta 全体にわたって導関数パターンを合成する。

    未知ビットがあっても, 全ての delta と全ての u で同じ値をとる出力ビットは
    確定する。ニブル単位で「非零なら全部未知」とするとこれを取りこぼし,
    フィルタ強度を過小評価する。
    """
    pats = [derivative_pattern(d) for d in deltas]
    out = []
    for k in range(8):
        vals = {p[k] for p in pats}
        out.append(vals.pop() if (len(vals) == 1 and UNK not in vals) else UNK)
    return out


def xor3(a, b):
    if a is UNK or b is UNK:
        return UNK
    return a ^ b


def f_layer(dx_arr, ncols, p_box):
    """
    ΔX から P(ΔG0 ⊕ ΔG1) を 3 値で計算する。

    入力ニブル m は G0 の出力ニブル m と G1 の出力ニブル (m-1) mod ncols を作る。
    ニブル m に未知ビットがあれば, そこから出る出力は全て未知とする。
    """
    nbits = 4 * ncols
    dg0 = [UNK] * nbits
    dg1 = [UNK] * nbits
    for m in range(ncols):
        pat = combined_pattern(possible_deltas(dx_arr, m, ncols))
        i = (m - 1) % ncols
        for k in range(4):                     # k=0 -> y3 (msb 行)
            dg0[k * ncols + m] = pat[k]
            dg1[k * ncols + i] = pat[4 + k]

    p_in = [xor3(dg0[t], dg1[t]) for t in range(nbits)]
    p_out = [UNK] * nbits
    for k in range(4):
        for i in range(ncols):
            p_out[i + ncols * k] = p_in[p_box[i] + ncols * k]
    return p_out


def truncated_backward(dx_s, dy_s, steps, ncols, p_box):
    """
    (ΔX^s, ΔY^s) が既知のとき, 復号方向に steps 段さかのぼった差分を 3 値で求める。

    復号方向のラウンド関係:
        ΔX^{r} = ΔY^{r+1}
        ΔY^{r} = ΔX^{r+1} ⊕ P(ΔG0(ΔX^r) ⊕ ΔG1(ΔX^r))

    返り値: [(dx_arr, dy_arr), ...]  index 0 が第 s 段, index t が第 s-t 段
    """
    nbits = 4 * ncols
    out = [(word_to_arr(dx_s, nbits), word_to_arr(dy_s, nbits))]
    for _ in range(steps):
        dx_next, dy_next = out[-1]
        dx_cur = list(dy_next)                                  # ΔX^r = ΔY^{r+1}
        p_out = f_layer(dx_cur, ncols, p_box)
        dy_cur = [xor3(dx_next[t], p_out[t]) for t in range(nbits)]
        out.append((dx_cur, dy_cur))
    return out


# --------------------------------------------------------------------------
def validate():
    """論文 (SAND-64) の Eq.(5), Eq.(6) を再現する"""
    ncols, p_box = 8, [7, 4, 1, 6, 3, 0, 5, 2]
    dx3, dy3 = 0x00900000, 0x008C0000
    seq = truncated_backward(dx3, dy3, 3, ncols, p_box)

    print("SAND-64, starting from (dX^3, dY^3) = (00900000, 008C0000), decryption direction")
    for t, (dxa, dya) in enumerate(seq):
        r = 3 - t
        print(f"  r={r}: dX = {arr_to_str(dxa, ncols)}   dY = {arr_to_str(dya, ncols)}")

    ok = True
    for r, expect_mask, expect_val, label in ((1, 0xB9FF37B9, 0x001E0000, "Eq.(6) dX^1"),
                                              (0, 0x81311081, 0x00110000, "Eq.(5) dX^0")):
        dxa, _ = seq[3 - r]
        m, v = known_mask(dxa, ncols)
        tag = "OK" if (m, v) == (expect_mask, expect_val) else "MISMATCH"
        if tag != "OK":
            ok = False
        print(f"  {label}: computed mask {m:08X} value {v:08X} (wt {bin(m).count('1')})"
              f"   paper {expect_mask:08X} / {expect_val:08X}  -> {tag}")
    print("validation PASSED" if ok else "validation FAILED")
    return ok


def analyse_sand128(dx_s, dy_s, steps=3):
    ncols, p_box = 16, K.P_BOX
    seq = truncated_backward(dx_s, dy_s, steps, ncols, p_box)
    print(f"SAND-128, from (dX^s, dY^s) = ({K.fmt(dx_s)}, {K.fmt(dy_s)}), "
          f"{steps} decryption rounds")
    filters = []
    for t, (dxa, dya) in enumerate(seq):
        m, v = known_mask(dxa, ncols)
        filters.append((t, m, v, bin(m).count("1")))
        print(f"  s-{t}: dX known bits {bin(m).count('1'):>3}  mask {m:016X}  value {v:016X}")
        print(f"        dX = {arr_to_str(dxa, ncols)}")
    return filters


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if cmd == "validate":
        validate()
    else:
        dx = int(sys.argv[2], 16)
        dy = int(sys.argv[3], 16)
        steps = int(sys.argv[4]) if len(sys.argv) > 4 else 3
        analyse_sand128(dx, dy, steps)
