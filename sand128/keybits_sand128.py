#!/usr/bin/env python3
"""
鍵回復に必要な丸鍵ビット数 k を厳密に数える (論文 Appendix A/B の一般化)。

差分では丸鍵が相殺するので, 鍵が要るのは「条件を検査するために状態の実値を
計算する」ところだけである。

    X^t = P(G0(X^{t-1}) ⊕ G1(X^{t-1})) ⊕ X^{t-2} ⊕ sk^{t-1}

なので X^t の座標 i を知るには sk^{t-1}[i] が要る。必要座標は 2 通りの経路で生じる:

  (1) 条件そのもの : ΔX^{c} を検査するには, ΔX^c = ΔX^{c-2} ⊕ P(D_{ΔX^{c-1}} f(X^{c-1}))
      より X^{c-1} の座標が要る。入力差分 delta のニブルについて必要なのは
      導関数が実際に依存する入力ビットだけ。ΔX^{c-1} に未知ビットがある場合は
      取り得る delta 全体にわたる和集合を取る (論文の M_K0^(1) がこれに当たる)。

  (2) 上位の必要座標の逆像 : X^t の座標 i を計算するには, P と G0/G1 を通じて
      i が依存する X^{t-1} の座標が要る (論文の M_K0^(2) がこれに当たる)。

検証: SAND-64 の論文の値
    sk^1 : 8C008C8C  (wt 9)
    sk^0 : DECE5EDE | D9D9D9D9 = DFDFDFDF  (wt 28)
    k = 37
"""
import sys

import filters_sand128 as F
import keyrecovery_sand128 as K


# --------------------------------------------------------------------------
# G0 / G1 の各出力ビットが依存する入力ビット (ANF から機械的に判定)
# --------------------------------------------------------------------------
def _func_dep(fn):
    """dep[b] = 出力ビット b が依存する入力ビットの集合"""
    dep = {}
    for b in range(4):
        s = set()
        for bp in range(4):
            if any(((fn(u) >> b) & 1) != ((fn(u ^ (1 << bp)) >> b) & 1) for u in range(16)):
                s.add(bp)
        dep[b] = s
    return dep


G = K.sand128.SAND_128_SAT_Model
G0DEP = _func_dep(G.g0_value)
G1DEP = _func_dep(G.g1_value)


def coord(ncols, b, col):
    """行 b (3 = msb), 列 col の座標を「index 0 = MSB」のビット位置に直す"""
    return (3 - b) * ncols + col


def mask_from_coords(coords, ncols):
    nbits = 4 * ncols
    m = 0
    for c in coords:
        m |= 1 << (nbits - 1 - c)
    return m


# --------------------------------------------------------------------------
# (1) 条件検査のために必要な状態座標
# --------------------------------------------------------------------------
def required_coords_from_diff(dx_arr, ncols):
    """
    差分 dx_arr (3 値配列) のニブルごとに, 導関数が依存する入力ビットを集める。
    未知ビットがある場合は取り得る delta 全体の和集合を取る。
    """
    out = set()
    for col in range(ncols):
        deltas = F.possible_deltas(dx_arr, col, ncols)
        sup = set()
        for d in deltas:
            sup |= K.derivative_support(d)
        for b in sup:
            out.add(coord(ncols, b, col))
    return out


# --------------------------------------------------------------------------
# (2) 上位の必要座標の逆像
# --------------------------------------------------------------------------
def preimage_coords(coords, ncols, p_box):
    """
    X^{t} の座標集合 coords を計算するために必要な X^{t-1} の座標集合。

    p_out[i + ncols*k] = p_in[p_box[i] + ncols*k] なので, 出力座標 (行 k, 列 i) は
    p_in の (行 k, 列 c=p_box[i]) に対応する。行 k は入出力ビット b = 3-k を表し,
    p_in の (b, c) は G0 の出力ニブル c のビット b と G1 の出力ニブル c のビット b の
    XOR である。G0 の出力ニブル c は入力ニブル c から, G1 の出力ニブル c は入力
    ニブル c+1 から作られる。
    """
    out = set()
    for c_idx in coords:
        k = c_idx // ncols
        i = c_idx % ncols
        b = 3 - k
        c = p_box[i]
        for bp in G0DEP[b]:
            out.add(coord(ncols, bp, c))
        for bp in G1DEP[b]:
            out.add(coord(ncols, bp, (c + 1) % ncols))
    return out


# --------------------------------------------------------------------------
# 全体
# --------------------------------------------------------------------------
def key_bits(dx_s, dy_s, r_in, ncols, p_box, verbose=False):
    """
    識別子を第 r_in 段に置き, 平文側に r_in 段付けたときに推測が必要な丸鍵ビット数。

    第 0, 1 段の差分条件は鍵なしで検査できる。第 2..r_in 段の条件検査に鍵が要る。
    返り値: (k, {t: sk^{t-1} のマスク})
    """
    seq = F.truncated_backward(dx_s, dy_s, r_in, ncols, p_box)
    # seq[t] は第 (r_in - t) 段。第 c 段の dX は seq[r_in - c]
    dx_of = {r_in - t: seq[t][0] for t in range(r_in + 1)}

    need = {t: set() for t in range(0, r_in + 1)}
    for c in range(2, r_in + 1):
        need[c - 1] |= required_coords_from_diff(dx_of[c - 1], ncols)
    for t in range(r_in - 1, 0, -1):
        need[t - 1] |= preimage_coords(need[t], ncols, p_box)

    masks = {}
    total = 0
    for t in range(1, r_in):
        if need[t]:
            m = mask_from_coords(need[t], ncols)
            masks[t] = m
            total += bin(m).count("1")
            if verbose:
                w = 4 * ncols // 4
                print(f"  sk^{t-1} mask = {m:0{ncols}X}  wt = {bin(m).count('1')}")
    return total, masks


def validate():
    """SAND-64 の論文の値 (sk^1 wt 9, sk^0 wt 28, k = 37) を再現する"""
    ncols, p_box = 8, [7, 4, 1, 6, 3, 0, 5, 2]
    dx3, dy3 = 0x00900000, 0x008C0000

    seq = F.truncated_backward(dx3, dy3, 3, ncols, p_box)
    dx2 = seq[1][0]
    dx1 = seq[2][0]

    m_k1 = mask_from_coords(required_coords_from_diff(dx2, ncols), ncols)
    print(f"sk^1 mask (from dX^2)                  = {m_k1:08X}  wt {bin(m_k1).count('1')}"
          f"   paper 8C008C8C / 9")
    ok = (m_k1 == 0x8C008C8C)

    m0a = mask_from_coords(required_coords_from_diff(dx1, ncols), ncols)
    print(f"sk^0 part 1 (from dX^1)                = {m0a:08X}  wt {bin(m0a).count('1')}"
          f"   paper DECE5EDE / 22")
    ok &= (m0a == 0xDECE5EDE)

    m0b = mask_from_coords(
        preimage_coords(required_coords_from_diff(dx2, ncols), ncols, p_box), ncols)
    print(f"sk^0 part 2 (preimage of sk^1 coords)  = {m0b:08X}  wt {bin(m0b).count('1')}"
          f"   paper D9D9D9D9")
    ok &= (m0b == 0xD9D9D9D9)

    m_k0 = m0a | m0b
    print(f"sk^0 mask (union)                      = {m_k0:08X}  wt {bin(m_k0).count('1')}"
          f"   paper DFDFDFDF / 28")
    ok &= (m_k0 == 0xDFDFDFDF)

    k, masks = key_bits(dx3, dy3, 3, ncols, p_box)
    print(f"total k = {k}   paper 37")
    ok &= (k == 37)
    print("validation PASSED" if ok else "validation FAILED")
    return ok


if __name__ == "__main__":
    if (sys.argv[1] if len(sys.argv) > 1 else "validate") == "validate":
        validate()
