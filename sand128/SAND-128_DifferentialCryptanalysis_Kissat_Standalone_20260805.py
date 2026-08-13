import subprocess
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.cnf import CnfMixin
from common.solver import SolverMixin
from common.ssb import SSbMixin
from common.tee import Tee

class SAND_128_SAT_Model(SSbMixin, CnfMixin, SolverMixin):
    def __init__(self):
        self.clauses = []
        self.vars = []
        self.p_box = [14,15,8,9,2,3,12,13,6,7,0,1,10,11,4,5]
        # Kissatの実行パス（リポジトリ直下の build/kissat が無ければ PATH 上の kissat を使う）
        _built = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "build", "kissat")
        self.kissat_path = _built if os.path.exists(_built) else (shutil.which("kissat") or _built)
        self.aux_var_counter = 0  # 補助変数のカウンター


    def generate_variables(self, num_rounds):
        variables = {}
        v = 0

        def next_var(name):
            nonlocal v
            v += 1
            variables[name] = v
            self.aux_var_counter = v  # 最大変数番号を記録

        # ラウンド0の入力
        for i in range(128):
            next_var(f"x_{i}_0")

        for round in range(1, num_rounds + 1):
            for i in range(128):
                next_var(f"x_{i}_{round}")
            for i in range(64):
                next_var(f"x_bottom_{i}_{round}")
            for i in range(64):
                next_var(f"g0_sx_{i}_{round}")
            for i in range(64):
                next_var(f"g0_sy_{i}_{round}")
            for i in range(64):
                next_var(f"g1_sx_{i}_{round}")
            for i in range(64):
                next_var(f"g1_sy_{i}_{round}")
            # 重みは G0/G1 ごとではなく, 両者を束ねた SSb ごとに 1 組 (unary 3 bit) 持つ
            for i in range(48):
                next_var(f"ssb_w_{i}_{round}")
            for i in range(64):
                next_var(f"p_input_{i}_{round}")

        return variables
        

    def s_box_model(self, variables, round):
        """
        SSb = (G0, G1) の同時差分制約を, 共有される入力ニブルごとに 1 個ずつ課す。

        ニブル j の入力差分は
          * G0 の j 番目の S-box    -> 出力ニブル j
          * G1 の (j-1) 番目の S-box -> 出力ニブル (j-1) mod 16
        の両方に入る。G1 側の入力が 1 ニブル左回転されているためである。
        G0 と G1 の出力差分は同じ入力値から生じなければならないので, 両者を
        単一の SSb 制約 (入力差分 4 bit + G0 出力差分 4 bit + G1 出力差分 4 bit
        + 重み) で結び付ける。以前のように周辺 DDT を独立に掛け合わせると,
        同時には実現し得ない遷移 (例: din = 0x6, (dG0, dG1) = (0x6, 0x2)) を
        許してしまう。

        g0_sx と g1_sx はいずれも branching により ΔX^r に等しいので, 入力差分は
        g0_sx から取れば足りる。
        """
        for j in range(16):
            a0, a1, a2, a3 = (variables[f"g0_sx_{j}_{round}"], variables[f"g0_sx_{j+16}_{round}"],
                              variables[f"g0_sx_{j+32}_{round}"], variables[f"g0_sx_{j+48}_{round}"])
            b0, b1, b2, b3 = (variables[f"g0_sy_{j}_{round}"], variables[f"g0_sy_{j+16}_{round}"],
                              variables[f"g0_sy_{j+32}_{round}"], variables[f"g0_sy_{j+48}_{round}"])

            # 入力ニブル j を消費する G1 の S-box は出力ニブル (j-1) mod 16 を作る
            i = (j - 1) % 16
            c0, c1, c2, c3 = (variables[f"g1_sy_{i}_{round}"], variables[f"g1_sy_{i+16}_{round}"],
                              variables[f"g1_sy_{i+32}_{round}"], variables[f"g1_sy_{i+48}_{round}"])

            w1, w2, w3 = (variables[f"ssb_w_{j}_{round}"], variables[f"ssb_w_{j+16}_{round}"],
                          variables[f"ssb_w_{j+32}_{round}"])

            for clause in self.ssb_s_box(a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, w1, w2, w3):
                self.clauses.append(clause)

            # Add the unary weight bits to self.vars for the sequential encoding method
            self.vars.extend([w1, w2, w3])

    def sequential_encoding_method(self, m):
        """
        Sequential counter encoding for at-most-m constraint
        This is a simple implementation of cardinality constraint
        """
        n = len(self.vars)
        if m >= n:
            return  # 制約なし
        
        # Sequential counter method
        # s[i][j]は、最初のi個の変数のうちj個以上が真であることを表す
        s = {}
        
        # 新しい補助変数を作成
        for i in range(1, n+1):
            for j in range(1, min(i, m+1)+1):
                self.aux_var_counter += 1
                s[(i,j)] = self.aux_var_counter
        
        # 制約を追加
        # s[1][1] = x[1]
        self.clauses.append([-self.vars[0], s[(1,1)]])
        self.clauses.append([self.vars[0], -s[(1,1)]])
        
        # 他の変数の初期化（範囲が正しくない場合はスキップ）
        for j in range(2, min(2, m+1)+1):
            if (1,j) in s:
                self.clauses.append([-s[(1,j)]])
        
        # 再帰的な定義
        for i in range(2, n+1):
            # s[i][1] = s[i-1][1] OR x[i]
            if (i-1,1) in s and (i,1) in s:
                self.clauses.append([-s[(i-1,1)], s[(i,1)]])
                self.clauses.append([-self.vars[i-1], s[(i,1)]])
                self.clauses.append([s[(i-1,1)], self.vars[i-1], -s[(i,1)]])
            
            for j in range(2, min(i, m+1)+1):
                # s[i][j] = s[i-1][j] OR (s[i-1][j-1] AND x[i])
                if (i-1,j) in s and (i,j) in s and (i-1,j-1) in s:
                    self.clauses.append([-s[(i-1,j)], s[(i,j)]])
                    self.clauses.append([-s[(i-1,j-1)], -self.vars[i-1], s[(i,j)]])
                    self.clauses.append([s[(i-1,j)], s[(i-1,j-1)], -s[(i,j)]])
                    self.clauses.append([s[(i-1,j)], self.vars[i-1], -s[(i,j)]])
                elif (i,j) in s and (i-1,j-1) in s:
                    # 対角要素 j = i では s[i-1][j] が存在しないので偽として扱う:
                    # s[i][i] = s[i-1][i-1] AND x[i]
                    # これを落とすと「最初の i 個がすべて真」という状況を検出できず,
                    # 制約が at-most-m ではなく at-most-(m+1) に緩んでしまう
                    self.clauses.append([-s[(i-1,j-1)], -self.vars[i-1], s[(i,j)]])
                    self.clauses.append([s[(i-1,j-1)], -s[(i,j)]])
                    self.clauses.append([self.vars[i-1], -s[(i,j)]])

        # 最後に、m+1個以上の変数が真にならないようにする
        if m+1 <= n and (n,m+1) in s:
            self.clauses.append([-s[(n,m+1)]])

    def hex_to_bits(self, hex_str):
        """64ビットの16進文字列をMSBファーストのビット列に変換する"""
        value = int(hex_str.replace(" ", ""), 16)
        return [(value >> (63 - i)) & 1 for i in range(64)]

    def bits_to_hex(self, bits):
        """MSBファーストの64ビット列を16進文字列に変換する"""
        value = 0
        for bit in bits:
            value = (value << 1) | bit
        return " ".join(f"{(value >> (48 - 16*k)) & 0xFFFF:04X}" for k in range(4))

    def fix_state_difference(self, variables, round, delta_x, delta_y):
        """
        ラウンドroundの状態差分 (ΔX^round, ΔY^round) を固定する
        ΔX^round: x_{i}_{round} (0 <= i < 64), ΔY^round: x_{i}_{round} (64 <= i < 128)
        変数の添字 i = 0 が MSB に対応する
        """
        for i, bit in enumerate(self.hex_to_bits(delta_x)):
            var = variables[f"x_{i}_{round}"]
            self.clauses.append([var] if bit == 1 else [-var])

        for i, bit in enumerate(self.hex_to_bits(delta_y)):
            var = variables[f"x_{i + 64}_{round}"]
            self.clauses.append([var] if bit == 1 else [-var])

    def permutation(self, input_bits, output_bits):
        """Permutation layer using the P-box."""
        for j in range(4):
            for i in range(16):
                output_bits[i+16*j] = input_bits[self.p_box[i]+16*j]


    def Search_Differential_Distingusher(self, num_rounds, num_bounds, variables):
        # Split the plaintext into two halves
        left_bits = [variables[str("x_" + str(i) + "_0")] for i in range(0,64)]
        right_bits = [variables[str("x_" + str(i) + "_0")] for i in range(64,128)]

        for round in range(1, num_rounds+1):
            # 1. Branching for G0
            g0_input = [variables[str("g0_sx_" + str(i) + "_" + str(round))] for i in range(0,64)]
            left_bottom = [variables[str("x_bottom_" + str(i) + "_" + str(round))] for i in range(0,64)]
            self.branching_model(left_bits, g0_input, left_bottom)

            # 2. Branching for G1
            g1_input = [variables[str("g1_sx_" + str(i) + "_" + str(round))] for i in range(0,64)]
            right_bottom = [variables[str("x_" + str(i) + "_" + str(round))] for i in range(64,128)]
            self.branching_model(left_bottom, g1_input, right_bottom)
            
            # 3. S-box transformation for G0 and G1
            self.s_box_model(variables, round)

            # 4. XOR operation between G0_output and G1_output
            g0_output = [variables[str("g0_sy_" + str(i) + "_" + str(round))] for i in range(0,64)]
            g1_output = [variables[str("g1_sy_" + str(i) + "_" + str(round))] for i in range(0,64)]
            p_input = [variables[str("p_input_" + str(i) + "_" + str(round))] for i in range(0,64)]
            self.xor_model(g0_output, g1_output, p_input)

            # 5. Permutation
            p_output = [0] * 64
            for j in range(4):
                for i in range(16):
                    p_output[i+16*j] = p_input[self.p_box[i]+16*j]

            # 5. XOR operation between p_output and right half
            right_after_xor = [variables[str("x_" + str(i) + "_" + str(round))] for i in range(0,64)]
            self.xor_model(p_output, right_bits, right_after_xor)

            # 6. Swap the halves for the next round
            left_bits = right_after_xor
            right_bits = right_bottom

        # Make sequential_encording model
        self.sequential_encoding_method(num_bounds)

    def Search_Differential_Distingusher_Backward(self, num_rounds, num_bounds, variables):
        """
        暗号文側から平文側への逆伝搬（復号方向）の差分特性を探索するモデル

        暗号化のラウンド関数は
            ΔX^{r+1} = ΔY^r xor P(G0(ΔX^r) xor G1(ΔX^r <<< 1)), ΔY^{r+1} = ΔX^r
        であるから, その逆（復号方向）は
            ΔX^r = ΔY^{r+1}, ΔY^r = ΔX^{r+1} xor P(G0(ΔX^r) xor G1(ΔX^r <<< 1))
        となる。

        逆伝搬のステップ数を t (1 <= t <= R) とし, 状態 (U^t, V^t) を
            (U^t, V^t) = (ΔX^{R-t}, ΔY^{R-t})
        と定義する。すなわち t = 0 が暗号文側 (ΔX^R, ΔY^R), t = R が平文側 (ΔX^0, ΔY^0)
        に対応する。変数 x_{i}_{t} は 0 <= i < 64 が U^t, 64 <= i < 128 が V^t を表す。
        """
        for t in range(1, num_rounds + 1):
            u_prev = [variables[str("x_" + str(i) + "_" + str(t-1))] for i in range(0,64)]
            v_prev = [variables[str("x_" + str(i) + "_" + str(t-1))] for i in range(64,128)]

            # 1. U^t = V^{t-1} であり, この値が G0, G1 の入力にも分岐する
            #    Branching for G0
            g0_input = [variables[str("g0_sx_" + str(i) + "_" + str(t))] for i in range(0,64)]
            bottom = [variables[str("x_bottom_" + str(i) + "_" + str(t))] for i in range(0,64)]
            self.branching_model(v_prev, g0_input, bottom)

            # 2. Branching for G1 and U^t
            g1_input = [variables[str("g1_sx_" + str(i) + "_" + str(t))] for i in range(0,64)]
            u_cur = [variables[str("x_" + str(i) + "_" + str(t))] for i in range(0,64)]
            self.branching_model(bottom, g1_input, u_cur)

            # 3. S-box transformation for G0 and G1
            self.s_box_model(variables, t)

            # 4. XOR operation between G0_output and G1_output
            g0_output = [variables[str("g0_sy_" + str(i) + "_" + str(t))] for i in range(0,64)]
            g1_output = [variables[str("g1_sy_" + str(i) + "_" + str(t))] for i in range(0,64)]
            p_input = [variables[str("p_input_" + str(i) + "_" + str(t))] for i in range(0,64)]
            self.xor_model(g0_output, g1_output, p_input)

            # 5. Permutation
            p_output = [0] * 64
            for j in range(4):
                for i in range(16):
                    p_output[i+16*j] = p_input[self.p_box[i]+16*j]

            # 6. V^t = XOR(U^{t-1}, P^t)
            v_cur = [variables[str("x_" + str(i) + "_" + str(t))] for i in range(64,128)]
            self.xor_model(p_output, u_prev, v_cur)

        # Make sequential_encording model
        self.sequential_encoding_method(num_bounds)


    def print_characteristics_v2(self, num_rounds, result_model, variables):
        def get_bit(var): return '1' if var in result_model else '0'

        # 入力ラウンドの状態ビット
        print("r = 0", end=" ")
        for i in range(128):
            if i == 64:
                print(end=" ")
            print(get_bit(variables[f"x_{i}_0"]), end="")
        print("")

        p=0
        # 各ラウンドの出力状態ビット
        for r in range(1, num_rounds + 1):
            print(f"r = {r}", end=" ")
            for i in range(128):
                if i == 64:
                    print(end=" ")
                print(get_bit(variables[f"x_{i}_{r}"]), end="")
            print(end=" ")

            # G0の状態ビット
            for i in range(64):
                print(get_bit(variables[f"g0_sx_{i}_{r}"]), end="")
            print(end=" ")
            for i in range(64):
                print(get_bit(variables[f"g0_sy_{i}_{r}"]), end="")
            print(end=" ")

            # G1の状態ビット
            for i in range(64):
                print(get_bit(variables[f"g1_sx_{i}_{r}"]), end="")
            print(end=" ")
            for i in range(64):
                print(get_bit(variables[f"g1_sy_{i}_{r}"]), end="")
            print(end=" ")

            # SSb ごとの重み (G0, G1 をまとめて 1 組)
            for i in range(48):
                print(get_bit(variables[f"ssb_w_{i}_{r}"]), end="")
                p += 1 * int(get_bit(variables[f"ssb_w_{i}_{r}"]))
            print(" p=2^-",p)

    def print_characteristics_forward(self, num_rounds, result_model, variables):
        """
        平文側 (r = 0) から暗号文側 (r = num_rounds) への差分特性を16進数で出力する

        変数 x_{i}_{r} は 0 <= i < 64 が ΔX^r, 64 <= i < 128 が ΔY^r を表す。
        log2(p) は平文側から数えた累積の重み（各ラウンドの SSb の重みの和）である。
        """
        model_set = set(result_model)
        def get_bit(var): return 1 if var in model_set else 0

        def state_hex(r):
            delta_x = [get_bit(variables[f"x_{i}_{r}"]) for i in range(0, 64)]
            delta_y = [get_bit(variables[f"x_{i}_{r}"]) for i in range(64, 128)]
            return self.bits_to_hex(delta_x), self.bits_to_hex(delta_y)

        print(f"{'r':>3} {'dX^r':>21} {'dY^r':>21} {'log2(p)':>9}   (r: round)")
        delta_x, delta_y = state_hex(0)
        print(f"{0:>3} {delta_x:>21} {delta_y:>21} {0:>9}   <- plaintext side")

        weight = 0
        for r in range(1, num_rounds + 1):
            # ラウンド r で通過する SSb の重み
            for i in range(48):
                weight += get_bit(variables[f"ssb_w_{i}_{r}"])
            delta_x, delta_y = state_hex(r)
            tail = "   <- ciphertext side" if r == num_rounds else ""
            print(f"{r:>3} {delta_x:>21} {delta_y:>21} {-weight:>9}{tail}")

        print(f"Differential probability of the {num_rounds}-round characteristic: p = 2^-{weight}")

    def print_characteristics_backward(self, num_rounds, result_model, variables):
        """
        暗号文側 (t = 0) から平文側 (t = num_rounds) への逆伝搬の差分特性を出力する
        t は逆伝搬のステップ数, r = num_rounds - t は暗号化のラウンド番号を表す
        """
        model_set = set(result_model)
        def get_bit(var): return 1 if var in model_set else 0

        def state_hex(t):
            delta_x = [get_bit(variables[f"x_{i}_{t}"]) for i in range(0, 64)]
            delta_y = [get_bit(variables[f"x_{i}_{t}"]) for i in range(64, 128)]
            return self.bits_to_hex(delta_x), self.bits_to_hex(delta_y)

        print(f"{'t':>3} {'r':>3} {'dX^r':>21} {'dY^r':>21} {'log2(p)':>9}   (t: backward step, r: round)")
        delta_x, delta_y = state_hex(0)
        print(f"{0:>3} {num_rounds:>3} {delta_x:>21} {delta_y:>21} {0:>9}   <- ciphertext side")

        weight = 0
        for t in range(1, num_rounds + 1):
            # 逆伝搬ステップ t で通過する SSb の重み（暗号化のラウンド num_rounds-t+1 に対応）
            for i in range(48):
                weight += get_bit(variables[f"ssb_w_{i}_{t}"])
            delta_x, delta_y = state_hex(t)
            tail = "   <- plaintext side" if t == num_rounds else ""
            print(f"{t:>3} {num_rounds - t:>3} {delta_x:>21} {delta_y:>21} {-weight:>9}{tail}")

        print(f"Differential probability of the {num_rounds}-round characteristic: p = 2^-{weight}")

    def run_backward_without_key(self, num_rounds, num_bounds, ciphertext_diff=None, plaintext_diff=None,
                                 max_plaintext_ones=None, max_output_ones=None):
        """
        選択暗号文攻撃の立場で, 暗号文側から平文側へ逆伝搬させる探索を行う
        ciphertext_diff  : (ΔX^R, ΔY^R)  … 逆伝搬の出発点 (t = 0) を固定する場合に指定する
        plaintext_diff   : (ΔX^0, ΔY^0)  … 平文側 (t = num_rounds) を固定する場合に指定する
        max_plaintext_ones: 平文差分 ΔY^0 の1の個数の上限（差分を固定しない場合の絞り込み用）
        max_output_ones   : 暗号文差分 ΔX^R の1の個数の上限（差分を固定しない場合の絞り込み用）
        """
        # define variables
        variables = self.generate_variables(num_rounds)

        # generate SAT model for the decryption (backward) direction
        self.Search_Differential_Distingusher_Backward(num_rounds, num_bounds, variables)

        # ciphertext differentials != all_0
        self.clauses.append([variables[str("x_" + str(i) + "_0")] for i in range(0,128)])

        # Fix the ciphertext difference (ΔX^R, ΔY^R) : starting point of the backward search
        if ciphertext_diff is not None:
            self.fix_state_difference(variables, 0, ciphertext_diff[0], ciphertext_diff[1])
            print(f"Fixed ciphertext difference: (dX^{num_rounds}, dY^{num_rounds}) = ({ciphertext_diff[0]}, {ciphertext_diff[1]})")

        # Fix the plaintext difference (ΔX^0, ΔY^0)
        if plaintext_diff is not None:
            self.fix_state_difference(variables, num_rounds, plaintext_diff[0], plaintext_diff[1])
            print(f"Fixed plaintext difference: (dX^0, dY^0) = ({plaintext_diff[0]}, {plaintext_diff[1]})")

        # 平文差分 ΔY^0 = x_{i}_{R} (64<=i<128) の1を最小化する
        if plaintext_diff is None and max_plaintext_ones is not None:
            self.vars = [variables[f"x_{i}_{num_rounds}"] for i in range(64, 128)]
            self.sequential_encoding_method(max_plaintext_ones)
            print(f"  - Plaintext differential dY^0 (x_i_{num_rounds} for 64<=i<128): max {max_plaintext_ones} ones")

        # 暗号文差分 ΔX^R = x_{i}_{0} (0<=i<64) の1を最小化する
        if ciphertext_diff is None and max_output_ones is not None:
            self.vars = [variables[f"x_{i}_0" ] for i in range(0, 64)]
            self.sequential_encoding_method(max_output_ones)
            print(f"  - Ciphertext differential dX^{num_rounds} (x_i_0 for 0<=i<64): max {max_output_ones} ones")

        # solve SAT model and print the result
        print(f"Solving SAND-128 differential cryptanalysis (backward) with {num_rounds} rounds and bound {num_bounds}")
        print("Using Kissat SAT solver (standalone version)...")
        print(f"Number of variables: {self.aux_var_counter}")
        print(f"Number of clauses: {len(self.clauses)}")

        result_boolean, result_model = self.solve()
        if result_boolean:
            print("The problem is satisfiable with model:")
            self.print_characteristics_backward(num_rounds, result_model, variables)
        else:
            print("The problem is unsatisfiable.")

    def run_without_key(self, num_rounds, num_bounds, max_plaintext_ones, max_output_ones):
        # define variables
        variables = self.generate_variables(num_rounds)

        # generate SAT model for SLIM without key
        self.Search_Differential_Distingusher(num_rounds, num_bounds, variables)
        
        # input differentials != all_0
        self.clauses.append([variables[str("x_" + str(i) + "_0")] for i in range(0,128)])

        # 平文差分の制約: x_{i}_0 (32<=i) の1を最小化
        # 右半分の64ビット（x_64_0からx_127_0）に対して制約を追加
        plaintext_vars_to_minimize = []
        for i in range(64, 128):  # 64から127まで（64ビット）
            plaintext_vars_to_minimize.append(variables[f"x_{i}_0"])
        
        # これらの変数の合計を制限
        self.vars = plaintext_vars_to_minimize
        self.sequential_encoding_method(max_plaintext_ones)
        
        # 最終段出力差分の制約: x_{i}_{r} (0<=i<64) の1を最小化
        output_vars_to_minimize = []
        for i in range(64):  # 0から63まで（64ビット）
            output_vars_to_minimize.append(variables[f"x_{i}_{num_rounds}"])
        
        # これらの変数の合計を制限
        self.vars = output_vars_to_minimize
        self.sequential_encoding_method(max_output_ones)

        # solve SAT model and print the result
        print(f"Solving SAND-128 differential cryptanalysis with {num_rounds} rounds and bound {num_bounds}")
        print("Using Kissat SAT solver (standalone version)...")
        print(f"Number of variables: {self.aux_var_counter}")
        print(f"Number of clauses: {len(self.clauses)}")
        print("Added constraints: ")
        print(f"  - Plaintext differential x_i_0 (64<=i<128): max {max_plaintext_ones} ones")
        print(f"  - Output differential x_i_{num_rounds} (0<=i<64): max {max_output_ones} ones")
        
        result_boolean, result_model = self.solve()
        if result_boolean:
            print("The problem is satisfiable with model:")
            self.print_characteristics_forward(num_rounds, result_model, variables)
        else:
            print("The problem is unsatisfiable.")

#--------------------------#
#       main function      #
#--------------------------#


if __name__ == "__main__":
    model = SAND_128_SAT_Model()
    # num_rounds:r, num_bounds:m
    r, b = 19, 120

    # 選択暗号文攻撃の立場: 暗号文側から平文側へ逆伝搬させる
    # 暗号文差分 (dX^19, dY^19) を固定する場合は次の行を設定する
    # ciphertext_diff = ("0000 0000 0000 0000", "0000 0000 0000 0000")
    ciphertext_diff = None
    # 平文差分 (dX^0, dY^0) を固定する場合は次の行を設定する
    # plaintext_diff = ("0000 0000 0000 0000", "0000 0000 0000 0000")
    plaintext_diff = None

    # 差分を固定しない側の絞り込み（1の個数の上限）
    max_plaintext_ones = 10   # 平文差分 dY^0 の1の個数の上限
    max_output_ones = 10      # 暗号文差分 dX^19 の1の個数の上限

    # 実行結果を results/results_SAND-128.txt に出力する
    _results = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(_results, exist_ok=True)
    result_file = os.path.join(_results, "results_SAND-128.txt")

    with open(result_file, 'w') as f:
        stdout = sys.stdout
        sys.stdout = Tee(stdout, f)
        try:
            print(f"Running SAND-128 differential cryptanalysis (forward / plaintext -> ciphertext):")
            print(f"Rounds: {r}, Bound: {b}")
            model.verify_ssb_model()
            print("-" * 50)
            model.run_without_key(r, b, max_plaintext_ones, max_output_ones)
            #model.run_backward_without_key(r, b, ciphertext_diff, plaintext_diff, max_plaintext_ones, max_output_ones)
        finally:
            sys.stdout = stdout

    print("-" * 50)
    print(f"Results are saved in {result_file}")

    # 順方向（平文側 -> 暗号文側）で探索する場合は run_without_key() を使用する
    # model.run_without_key(r, b, max_plaintext_ones, max_output_ones)
