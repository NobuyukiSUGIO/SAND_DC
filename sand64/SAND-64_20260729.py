import subprocess
import os
import shutil
import sys
import tempfile

class SAND_64_SAT_Model():
    def __init__(self):
        self.clauses = []
        self.vars = []
        self.p_box = [7,4,1,6,3,0,5,2]
        # Kissatの実行パス（リポジトリ直下の build/kissat が無ければ PATH 上の kissat を使う）
        _built = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "build", "kissat")
        self.kissat_path = _built if os.path.exists(_built) else (shutil.which("kissat") or _built)
        self.aux_var_counter = 0  # 補助変数のカウンター

    @staticmethod
    def g0_value(x):
        """G0 acting on a 4-bit value x (bit 3 = msb). Used only for verification."""
        x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
        y0 = (x3 & x2) ^ x0
        y1 = x1
        y2 = x2
        y3 = (y0 & x1) ^ x3
        return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0

    @staticmethod
    def g1_value(x):
        """G1 acting on a 4-bit value x (bit 3 = msb). Used only for verification."""
        x3, x2, x1, x0 = (x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1
        y3 = x3
        y2 = (x3 & x1) ^ x2
        y1 = (y2 & x0) ^ x1
        y0 = x0
        return (y3 << 3) | (y2 << 2) | (y1 << 1) | y0

    @classmethod
    def ssb_joint_ddt(cls):
        """
        Differential distribution table of the 4x8 synthetic S-box

            SSb(x) = G0(x) || G1(x),

        i.e. ddt[din][(dG0 << 4) | dG1] counts the 4-bit values x with
        G0(x) ^ G0(x ^ din) = dG0 and G1(x) ^ G1(x ^ din) = dG1 simultaneously.

        G0 and G1 consume the *same* 4-bit input value, so their output
        differences must not be drawn from two independent DDTs: the product of
        the marginal DDTs admits 157 (din, dG0, dG1) triples whereas only 97 are
        jointly realizable.
        """
        ddt = [[0] * 256 for _ in range(16)]
        for din in range(16):
            for x in range(16):
                d0 = cls.g0_value(x) ^ cls.g0_value(x ^ din)
                d1 = cls.g1_value(x) ^ cls.g1_value(x ^ din)
                ddt[din][(d0 << 4) | d1] += 1
        return ddt

    def ssb_s_box(self, a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, w1, w2, w3):
        """
        Joint differential constraint for the synthetic S-box SSb = (G0, G1).

        Number of constraints: 58
        Input:      a0||a1||a2||a3; a0: msb   (input difference of the shared nibble)
        Output G0:  b0||b1||b2||b3; b0: msb
        Output G1:  c0||c1||c2||c3; c0: msb
        Weight:     w1 + w2 + w3, encoded in unary with w1 >= w2 >= w3.

        The nonzero entries of the joint DDT are 2, 4 and 16, so the weight
        -log2(Pr) of a valid transition is 3, 2 or 0; a weight of 1 never occurs
        and is excluded by the clauses below.  The clause set was generated from
        ssb_joint_ddt() and its models were exhaustively checked against the
        table by verify_ssb_model().
        """
        return [[-a0,c0],[-a0,w2],[-a1,b1],[-a2,b2],[-a3,c3],[-b3,w2],[-c1,w1],[-c2,w1],[-w1,w2],[a0,-c0],[a1,-b1],[a2,-b2],[a3,-c3],[w1,-w2],[w1,-w3],[-c0,-c1,w3],[b2,c1,-w3],[-a0,b0,b2,b3],[-a2,-a3,c1,w3],[-a2,c1,c2,c3],[-b0,-c1,-c2,w3],[-b0,b2,b3,c0],[-b1,-b2,-c3,w3],[a1,b0,-c0,w3],[a2,-b1,c0,c1],[a3,-c0,-c2,w3],[b0,-b2,c2,w3],[b1,-b3,c0,c3],[b1,b3,-c1,w3],[b1,b3,c0,-c3],[b2,-b3,-c1,w3],[c0,c1,c3,-w3],[-a0,a2,b0,-b1,-c1],[-a1,-a3,-b3,c2,w3],[-a1,-b3,c0,c1,-w3],[-a1,b3,c0,-c1,-w3],[-a2,-b3,c0,-c1,c3],[-a3,b1,c0,-c1,-w3],[-b0,-b3,-c0,-c2,w3],[-b0,b1,b2,-b3,-c1],[-b1,b3,c1,-c2,w3],[a1,-a3,b3,c2,w3],[a1,b0,b2,b3,-w1],[a2,-b0,-b1,-b3,w3],[-a0,-b0,-b2,b3,c2,c3],[-a0,b0,b1,b3,-c1,-c2],[-a1,-a2,b3,-c0,c1,c2],[-a1,-b2,-b3,c1,-c2,-c3],[-a1,b0,-c0,-c1,-c2,c3],[-a3,b1,-b2,b3,c1,-c2],[-b0,b1,-b3,-c1,-c2,c3],[a1,-b3,-c0,c1,c2,-w3],[a1,a3,b0,-b2,-b3,c2],[b0,-b1,b3,-c0,c2,-c3],[-a0,-a1,-a2,-b0,-b3,-c1,c2],[-a0,b0,-b2,-b3,-c1,-c2,-c3],[-a1,-a2,-b0,b3,-c1,-c2,-c3],[-a2,-b0,b1,-c0,-c1,c2,-c3]]

    def verify_ssb_model(self, verbose=True):
        """
        Exhaustively check that the models of ssb_s_box() are exactly the valid
        (input difference, G0 output difference, G1 output difference, weight)
        tuples of the joint DDT.  All 2^15 assignments of the 15 involved
        variables are enumerated, so this is a proof, not a sample.
        """
        clauses = self.ssb_s_box(*range(1, 16))  # variables 1..15 in signature order
        ddt = self.ssb_joint_ddt()

        expected = set()
        for din in range(16):
            for dout, count in enumerate(ddt[din]):
                if count == 0:
                    continue
                weight = 4 - count.bit_length() + 1
                assert count == 1 << (count.bit_length() - 1), "non power-of-two DDT entry"
                bits = []
                for value in (din, (dout >> 4) & 0xF, dout & 0xF):
                    bits += [(value >> 3) & 1, (value >> 2) & 1, (value >> 1) & 1, value & 1]
                bits += [1 if k < weight else 0 for k in range(3)]
                expected.add(tuple(bits))

        models = set()
        for code in range(1 << 15):
            bits = [(code >> (14 - k)) & 1 for k in range(15)]
            assign = {k + 1: bits[k] for k in range(15)}
            if all(any(assign[abs(l)] == (1 if l > 0 else 0) for l in clause) for clause in clauses):
                models.add(tuple(bits))

        if models != expected:
            raise AssertionError(
                f"SSb CNF is not exact: {len(models - expected)} spurious and "
                f"{len(expected - models)} missing assignments"
            )
        if verbose:
            print(f"SSb joint-DDT model verified: {len(expected)} valid transitions, "
                  f"{len(clauses)} clauses, all 2^15 assignments checked")
        return True


    def generate_variables(self, num_rounds):
        variables = {}
        v = 0

        def next_var(name):
            nonlocal v
            v += 1
            variables[name] = v
            self.aux_var_counter = v  # 最大変数番号を記録

        # ラウンド0の入力
        for i in range(64):
            next_var(f"x_{i}_0")

        for round in range(1, num_rounds + 1):
            for i in range(64):
                next_var(f"x_{i}_{round}")
            for i in range(32):
                next_var(f"x_bottom_{i}_{round}")
            for i in range(32):
                next_var(f"g0_sx_{i}_{round}")
            for i in range(32):
                next_var(f"g0_sy_{i}_{round}")
            for i in range(32):
                next_var(f"g1_sx_{i}_{round}")
            for i in range(32):
                next_var(f"g1_sy_{i}_{round}")
            # 重みは G0/G1 ごとではなく, 両者を束ねた SSb ごとに 1 組 (unary 3 bit) 持つ
            for i in range(24):
                next_var(f"ssb_w_{i}_{round}")
            for i in range(32):
                next_var(f"p_input_{i}_{round}")

        return variables
        
    def rotate_left(self, data, bits):
        """Rotate a list to the left by a given number of bits."""
        return data[bits:] + data[:bits]

    def s_box_model(self, variables, round):
        """
        SSb = (G0, G1) の同時差分制約を, 共有される入力ニブルごとに 1 個ずつ課す。

        ニブル j の入力差分は
          * G0 の j 番目の S-box   -> 出力ニブル j
          * G1 の (j-1) 番目の S-box -> 出力ニブル (j-1) mod 8
        の両方に入る。G1 側の入力が 1 ニブル左回転されているためである。
        G0 と G1 の出力差分は同じ入力値から生じなければならないので, 両者を
        単一の SSb 制約 (入力差分 4 bit + G0 出力差分 4 bit + G1 出力差分 4 bit
        + 重み) で結び付ける。以前のように周辺 DDT を独立に掛け合わせると,
        同時には実現し得ない遷移 (例: din = 0x6, (dG0, dG1) = (0x6, 0x2)) を
        許してしまう。

        g0_sx と g1_sx はいずれも branching により ΔX^r に等しいので, 入力差分は
        g0_sx から取れば足りる。
        """
        for j in range(8):
            a0, a1, a2, a3 = (variables[f"g0_sx_{j}_{round}"], variables[f"g0_sx_{j+8}_{round}"],
                              variables[f"g0_sx_{j+16}_{round}"], variables[f"g0_sx_{j+24}_{round}"])
            b0, b1, b2, b3 = (variables[f"g0_sy_{j}_{round}"], variables[f"g0_sy_{j+8}_{round}"],
                              variables[f"g0_sy_{j+16}_{round}"], variables[f"g0_sy_{j+24}_{round}"])

            # 入力ニブル j を消費する G1 の S-box は出力ニブル (j-1) mod 8 を作る
            i = (j - 1) % 8
            c0, c1, c2, c3 = (variables[f"g1_sy_{i}_{round}"], variables[f"g1_sy_{i+8}_{round}"],
                              variables[f"g1_sy_{i+16}_{round}"], variables[f"g1_sy_{i+24}_{round}"])

            w1, w2, w3 = (variables[f"ssb_w_{j}_{round}"], variables[f"ssb_w_{j+8}_{round}"],
                          variables[f"ssb_w_{j+16}_{round}"])

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
    
    def sequential_encoding_method_for_vars(self, vars_list, m):
        """
        Sequential counter encoding for at-most-m constraint on specific variables
        """
        n = len(vars_list)
        if m >= n:
            return  # 制約なし
        
        # Sequential counter method
        s = {}
        
        # 新しい補助変数を作成
        for i in range(1, n+1):
            for j in range(1, min(i, m+1)+1):
                self.aux_var_counter += 1
                s[(i,j)] = self.aux_var_counter
        
        # 制約を追加
        # s[1][1] = vars_list[0]
        self.clauses.append([-vars_list[0], s[(1,1)]])
        self.clauses.append([vars_list[0], -s[(1,1)]])
        
        # 他の変数の初期化
        for j in range(2, min(2, m+1)+1):
            if (1,j) in s:
                self.clauses.append([-s[(1,j)]])
        
        # 再帰的な定義
        for i in range(2, n+1):
            # s[i][1] = s[i-1][1] OR vars_list[i-1]
            if (i-1,1) in s and (i,1) in s:
                self.clauses.append([-s[(i-1,1)], s[(i,1)]])
                self.clauses.append([-vars_list[i-1], s[(i,1)]])
                self.clauses.append([s[(i-1,1)], vars_list[i-1], -s[(i,1)]])
            
            for j in range(2, min(i, m+1)+1):
                # s[i][j] = s[i-1][j] OR (s[i-1][j-1] AND vars_list[i-1])
                if (i-1,j) in s and (i,j) in s and (i-1,j-1) in s:
                    self.clauses.append([-s[(i-1,j)], s[(i,j)]])
                    self.clauses.append([-s[(i-1,j-1)], -vars_list[i-1], s[(i,j)]])
                    self.clauses.append([s[(i-1,j)], s[(i-1,j-1)], -s[(i,j)]])
                    self.clauses.append([s[(i-1,j)], vars_list[i-1], -s[(i,j)]])
                elif (i,j) in s and (i-1,j-1) in s:
                    # 対角要素 j = i: s[i][i] = s[i-1][i-1] AND vars_list[i-1]
                    self.clauses.append([-s[(i-1,j-1)], -vars_list[i-1], s[(i,j)]])
                    self.clauses.append([s[(i-1,j-1)], -s[(i,j)]])
                    self.clauses.append([vars_list[i-1], -s[(i,j)]])

        # 最後に、m+1個以上の変数が真にならないようにする
        if m+1 <= n and (n,m+1) in s:
            self.clauses.append([-s[(n,m+1)]])

    def hex_to_bits(self, hex_str):
        """32ビットの16進文字列をMSBファーストのビット列に変換する"""
        value = int(hex_str.replace(" ", ""), 16)
        return [(value >> (31 - i)) & 1 for i in range(32)]

    def bits_to_hex(self, bits):
        """MSBファーストの32ビット列を16進文字列に変換する"""
        value = 0
        for bit in bits:
            value = (value << 1) | bit
        return f"{value >> 16:04X} {value & 0xFFFF:04X}"

    def fix_state_difference(self, variables, round, delta_x, delta_y):
        """
        ラウンドroundの状態差分 (ΔX^round, ΔY^round) を固定する
        ΔX^round: x_{i}_{round} (0 <= i < 32), ΔY^round: x_{i}_{round} (32 <= i < 64)
        変数の添字 i = 0 が MSB に対応する
        """
        for i, bit in enumerate(self.hex_to_bits(delta_x)):
            var = variables[f"x_{i}_{round}"]
            self.clauses.append([var] if bit == 1 else [-var])

        for i, bit in enumerate(self.hex_to_bits(delta_y)):
            var = variables[f"x_{i + 32}_{round}"]
            self.clauses.append([var] if bit == 1 else [-var])

    def permutation(self, input_bits, output_bits):
        """Permutation layer using the P-box."""
        for j in range(4):
            for i in range(8):
                output_bits[i+8*j] = input_bits[self.p_box[i]+8*j]

    def xor_model(self, input1_bits, input2_bits, output_bits):
        # For each bit position, model the XOR operation
        for alpha, beta, gamma in zip(input1_bits, input2_bits, output_bits):
            # If alpha and beta are both 0, then gamma should be 0
            self.clauses.append([-alpha, -beta, -gamma])
            # If alpha is 1 and beta is 0, then gamma should be 1
            self.clauses.append([alpha, -beta, gamma])
            # If alpha is 0 and beta is 1, then gamma should be 1
            self.clauses.append([-alpha, beta, gamma])
            # If alpha and beta are both 1, then gamma should be 0
            self.clauses.append([alpha, beta, -gamma])

    def branching_model(self, input_bits, output1_bits, output2_bits):
        # input_bits: 入力ビットのリスト
        # output1_bits: 出力ビットのリスト1
        # output2_bits: 出力ビットのリスト2

        # ブランチング操作のモデル化
        for i in range(len(input_bits)):
            # 入力ビットが1の場合、出力ビット1と出力ビット2の両方が1である必要がある
            self.clauses.append([-input_bits[i], output1_bits[i]])
            self.clauses.append([-input_bits[i], output2_bits[i]])

            # 入力ビットが0の場合、出力ビット1と出力ビット2の両方が0である必要がある
            self.clauses.append([input_bits[i], -output1_bits[i]])
            self.clauses.append([input_bits[i], -output2_bits[i]])

    def copy_model(self, input, output):
        # input_bits: 入力ビット
        # output_bits: 出力ビット

        # 入力ビットが1の場合、出力ビットが1である必要がある
        self.clauses.append([-input, output])

        # 入力ビットが0の場合、出力ビットが0である必要がある
        self.clauses.append([input, -output])

    def Search_Differential_Distingusher(self, num_rounds, num_bounds, variables):
        # Split the plaintext into two halves
        left_bits = [variables[str("x_" + str(i) + "_0")] for i in range(0,32)]
        right_bits = [variables[str("x_" + str(i) + "_0")] for i in range(32,64)]

        for round in range(1, num_rounds+1):
            # 1. Branching for G0
            g0_input = [variables[str("g0_sx_" + str(i) + "_" + str(round))] for i in range(0,32)]
            left_bottom = [variables[str("x_bottom_" + str(i) + "_" + str(round))] for i in range(0,32)]
            self.branching_model(left_bits, g0_input, left_bottom)

            # 2. Branching for G1
            g1_input = [variables[str("g1_sx_" + str(i) + "_" + str(round))] for i in range(0,32)]
            right_bottom = [variables[str("x_" + str(i) + "_" + str(round))] for i in range(32,64)]
            self.branching_model(left_bottom, g1_input, right_bottom)
            
            # 3. S-box transformation for G0 and G1
            self.s_box_model(variables, round)

            # 4. XOR operation between G0_output and G1_output
            g0_output = [variables[str("g0_sy_" + str(i) + "_" + str(round))] for i in range(0,32)]
            g1_output = [variables[str("g1_sy_" + str(i) + "_" + str(round))] for i in range(0,32)]
            p_input = [variables[str("p_input_" + str(i) + "_" + str(round))] for i in range(0,32)]
            self.xor_model(g0_output, g1_output, p_input)

            # 5. Permutation
            p_output = [0] * 32
            for j in range(4):
                for i in range(8):
                    p_output[i+8*j] = p_input[self.p_box[i]+8*j]

            # 5. XOR operation between p_output and right half
            right_after_xor = [variables[str("x_" + str(i) + "_" + str(round))] for i in range(0,32)]
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
        に対応する。変数 x_{i}_{t} は 0 <= i < 32 が U^t, 32 <= i < 64 が V^t を表す。
        """
        for t in range(1, num_rounds + 1):
            u_prev = [variables[str("x_" + str(i) + "_" + str(t-1))] for i in range(0,32)]
            v_prev = [variables[str("x_" + str(i) + "_" + str(t-1))] for i in range(32,64)]

            # 1. U^t = V^{t-1} であり, この値が G0, G1 の入力にも分岐する
            #    Branching for G0
            g0_input = [variables[str("g0_sx_" + str(i) + "_" + str(t))] for i in range(0,32)]
            bottom = [variables[str("x_bottom_" + str(i) + "_" + str(t))] for i in range(0,32)]
            self.branching_model(v_prev, g0_input, bottom)

            # 2. Branching for G1 and U^t
            g1_input = [variables[str("g1_sx_" + str(i) + "_" + str(t))] for i in range(0,32)]
            u_cur = [variables[str("x_" + str(i) + "_" + str(t))] for i in range(0,32)]
            self.branching_model(bottom, g1_input, u_cur)

            # 3. S-box transformation for G0 and G1
            self.s_box_model(variables, t)

            # 4. XOR operation between G0_output and G1_output
            g0_output = [variables[str("g0_sy_" + str(i) + "_" + str(t))] for i in range(0,32)]
            g1_output = [variables[str("g1_sy_" + str(i) + "_" + str(t))] for i in range(0,32)]
            p_input = [variables[str("p_input_" + str(i) + "_" + str(t))] for i in range(0,32)]
            self.xor_model(g0_output, g1_output, p_input)

            # 5. Permutation
            p_output = [0] * 32
            for j in range(4):
                for i in range(8):
                    p_output[i+8*j] = p_input[self.p_box[i]+8*j]

            # 6. V^t = XOR(U^{t-1}, P^t)
            v_cur = [variables[str("x_" + str(i) + "_" + str(t))] for i in range(32,64)]
            self.xor_model(p_output, u_prev, v_cur)

        # Make sequential_encording model
        self.sequential_encoding_method(num_bounds)

    def write_cnf_to_file(self, filename):
        """CNF式をDIMACS形式でファイルに書き込む"""
        with open(filename, 'w') as f:
            # ヘッダーを書き込む
            num_vars = max(max(abs(lit) for lit in clause) for clause in self.clauses)
            num_clauses = len(self.clauses)
            f.write(f"p cnf {num_vars} {num_clauses}\n")
            
            # 各節を書き込む
            for clause in self.clauses:
                f.write(' '.join(map(str, clause)) + ' 0\n')

    def solve_with_kissat(self, cnf_file):
        """Kissatを使用してCNFを解く（時間制限なし）"""
        try:
            # Kissatを実行（時間制限を排除）
            result = subprocess.run([self.kissat_path, cnf_file], 
                                    capture_output=True, text=True)
            
            # 結果を解析
            output = result.stdout
            if "s SATISFIABLE" in output:
                # モデルを抽出
                model = []
                for line in output.split('\n'):
                    if line.startswith('v '):
                        vars_str = line[2:].strip()
                        for var in vars_str.split():
                            if var == '0':
                                break
                            model.append(int(var))
                return True, model
            elif "s UNSATISFIABLE" in output:
                return False, []
            else:
                print("Kissat output:")
                print(output)
                return False, []
                
        except Exception as e:
            print(f"Error running Kissat: {e}")
            return False, []

    def solve(self):
        # 一時ファイルを作成してCNFを書き込む
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cnf', delete=False) as f:
            temp_filename = f.name
        
        try:
            # CNFをファイルに書き込む
            self.write_cnf_to_file(temp_filename)
            
            # Kissatで解く
            sat, model = self.solve_with_kissat(temp_filename)
            
            # モデルを正の値のみのリストに変換（pysat形式に合わせる）
            if sat:
                result_model = [var for var in model if var > 0]
            else:
                result_model = []
                
            return sat, result_model
            
        finally:
            # 一時ファイルを削除
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    def print_characteristics_v2(self, num_rounds, result_model, variables):
        def get_bit(var): return '1' if var in result_model else '0'

        # 入力ラウンドの状態ビット
        print("r = 0", end=" ")
        for i in range(64):
            if i == 32:
                print(end=" ")
            print(get_bit(variables[f"x_{i}_0"]), end="")
        print("")

        p=0
        # 各ラウンドの出力状態ビット
        for r in range(1, num_rounds + 1):
            print(f"r = {r}", end=" ")
            for i in range(64):
                if i == 32:
                    print(end=" ")
                print(get_bit(variables[f"x_{i}_{r}"]), end="")
            print(end=" ")

            # G0の状態ビット
            for i in range(32):
                print(get_bit(variables[f"g0_sx_{i}_{r}"]), end="")
            print(end=" ")
            for i in range(32):
                print(get_bit(variables[f"g0_sy_{i}_{r}"]), end="")
            print(end=" ")

            # G1の状態ビット
            for i in range(32):
                print(get_bit(variables[f"g1_sx_{i}_{r}"]), end="")
            print(end=" ")
            for i in range(32):
                print(get_bit(variables[f"g1_sy_{i}_{r}"]), end="")
            print(end=" ")

            # SSb ごとの重み (G0, G1 をまとめて 1 組)
            for i in range(24):
                print(get_bit(variables[f"ssb_w_{i}_{r}"]), end="")
                p += 1 * int(get_bit(variables[f"ssb_w_{i}_{r}"]))
            print(" p=2^-",p)

    def print_characteristics_forward(self, num_rounds, result_model, variables):
        """
        平文側 (r = 0) から暗号文側 (r = num_rounds) への差分特性を16進数で出力する

        変数 x_{i}_{r} は 0 <= i < 32 が ΔX^r, 32 <= i < 64 が ΔY^r を表す。
        log2(p) は平文側から数えた累積の重み（各ラウンドの SSb の重みの和）である。
        """
        model_set = set(result_model)
        def get_bit(var): return 1 if var in model_set else 0

        def state_hex(r):
            delta_x = [get_bit(variables[f"x_{i}_{r}"]) for i in range(0, 32)]
            delta_y = [get_bit(variables[f"x_{i}_{r}"]) for i in range(32, 64)]
            return self.bits_to_hex(delta_x), self.bits_to_hex(delta_y)

        print(f"{'r':>3} {'dX^r':>11} {'dY^r':>11} {'log2(p)':>9}   (r: round)")
        delta_x, delta_y = state_hex(0)
        print(f"{0:>3} {delta_x:>11} {delta_y:>11} {0:>9}   <- plaintext side")

        weight = 0
        for r in range(1, num_rounds + 1):
            # ラウンド r で通過する SSb の重み
            for i in range(24):
                weight += get_bit(variables[f"ssb_w_{i}_{r}"])
            delta_x, delta_y = state_hex(r)
            tail = "   <- ciphertext side" if r == num_rounds else ""
            print(f"{r:>3} {delta_x:>11} {delta_y:>11} {-weight:>9}{tail}")

        print(f"Differential probability of the {num_rounds}-round characteristic: p = 2^-{weight}")

    def print_characteristics_backward(self, num_rounds, result_model, variables):
        """
        暗号文側 (t = 0) から平文側 (t = num_rounds) への逆伝搬の差分特性を出力する
        t は逆伝搬のステップ数, r = num_rounds - t は暗号化のラウンド番号を表す
        """
        model_set = set(result_model)
        def get_bit(var): return 1 if var in model_set else 0

        def state_hex(t):
            delta_x = [get_bit(variables[f"x_{i}_{t}"]) for i in range(0, 32)]
            delta_y = [get_bit(variables[f"x_{i}_{t}"]) for i in range(32, 64)]
            return self.bits_to_hex(delta_x), self.bits_to_hex(delta_y)

        print(f"{'t':>3} {'r':>3} {'dX^r':>11} {'dY^r':>11} {'log2(p)':>9}   (t: backward step, r: round)")
        delta_x, delta_y = state_hex(0)
        print(f"{0:>3} {num_rounds:>3} {delta_x:>11} {delta_y:>11} {0:>9}   <- ciphertext side (fixed)")

        weight = 0
        for t in range(1, num_rounds + 1):
            # 逆伝搬ステップ t で通過する SSb の重み（暗号化のラウンド num_rounds-t+1 に対応）
            for i in range(24):
                weight += get_bit(variables[f"ssb_w_{i}_{t}"])
            delta_x, delta_y = state_hex(t)
            tail = "   <- plaintext side" if t == num_rounds else ""
            print(f"{t:>3} {num_rounds - t:>3} {delta_x:>11} {delta_y:>11} {-weight:>9}{tail}")

        print(f"Differential probability of the {num_rounds}-round characteristic: p = 2^-{weight}")

    def run_backward_without_key(self, num_rounds, num_bounds, ciphertext_diff=None, plaintext_diff=None):
        """
        選択暗号文攻撃の立場で, 暗号文側の差分を固定して平文側へ逆伝搬させる探索を行う
        ciphertext_diff: (ΔX^R, ΔY^R)  … 逆伝搬の出発点 (t = 0)
        plaintext_diff : (ΔX^0, ΔY^0)  … 平文側 (t = num_rounds) も固定する場合に指定する
        """
        # define variables
        variables = self.generate_variables(num_rounds)

        # generate SAT model for the decryption (backward) direction
        self.Search_Differential_Distingusher_Backward(num_rounds, num_bounds, variables)

        # ciphertext differentials != all_0
        self.clauses.append([variables[str("x_" + str(i) + "_0")] for i in range(0,64)])

        # Fix the ciphertext difference (ΔX^R, ΔY^R) : starting point of the backward search
        if ciphertext_diff is not None:
            self.fix_state_difference(variables, 0, ciphertext_diff[0], ciphertext_diff[1])
            print(f"Fixed ciphertext difference: (dX^{num_rounds}, dY^{num_rounds}) = ({ciphertext_diff[0]}, {ciphertext_diff[1]})")

        # Fix the plaintext difference (ΔX^0, ΔY^0)
        if plaintext_diff is not None:
            self.fix_state_difference(variables, num_rounds, plaintext_diff[0], plaintext_diff[1])
            print(f"Fixed plaintext difference: (dX^0, dY^0) = ({plaintext_diff[0]}, {plaintext_diff[1]})")

        # solve SAT model and print the result
        print(f"Solving SAND-64 differential cryptanalysis (backward) with {num_rounds} rounds and bound {num_bounds}")
        print("Using Kissat SAT solver (standalone version)...")
        print(f"Number of variables: {self.aux_var_counter}")
        print(f"Number of clauses: {len(self.clauses)}")

        result_boolean, result_model = self.solve()
        if result_boolean:
            print("The problem is satisfiable with model:")
            self.print_characteristics_backward(num_rounds, result_model, variables)
        else:
            print("The problem is unsatisfiable.")

    def run_without_key(self, num_rounds, num_bounds, output_diff=None, input_diff=None):
        """
        平文側から暗号文側へ順方向に伝搬させる探索を行う
        output_diff: (ΔX^R, ΔY^R) … 暗号文側 (r = num_rounds) を固定する場合に指定する
        input_diff : (ΔX^0, ΔY^0) … 平文側 (r = 0) を固定する場合に指定する
        差分特性は print_characteristics_forward() により16進数で出力される
        """
        # define variables
        variables = self.generate_variables(num_rounds)

        # generate SAT model for SLIM without key
        self.Search_Differential_Distingusher(num_rounds, num_bounds, variables)

        # input differentials != all_0
        self.clauses.append([variables[str("x_" + str(i) + "_0")] for i in range(0,64)])

        # Fix the input difference (ΔX^0, ΔY^0)
        if input_diff is not None:
            self.fix_state_difference(variables, 0, input_diff[0], input_diff[1])
            print(f"Fixed input difference: (dX^0, dY^0) = ({input_diff[0]}, {input_diff[1]})")

        # Fix the output difference of the last round (ΔX^r, ΔY^r)
        if output_diff is not None:
            self.fix_state_difference(variables, num_rounds, output_diff[0], output_diff[1])
            print(f"Fixed output difference: (dX^{num_rounds}, dY^{num_rounds}) = ({output_diff[0]}, {output_diff[1]})")

        # solve SAT model and print the result
        print(f"Solving SAND-64 differential cryptanalysis with {num_rounds} rounds and bound {num_bounds}")
        print("Using Kissat SAT solver (standalone version)...")
        print(f"Number of variables: {self.aux_var_counter}")
        print(f"Number of clauses: {len(self.clauses)}")
        
        result_boolean, result_model = self.solve()
        if result_boolean:
            print("The problem is satisfiable with model:")
            self.print_characteristics_forward(num_rounds, result_model, variables)
        else:
            print("The problem is unsatisfiable.")

#--------------------------#
#       main function      #
#--------------------------#

class Tee:
    """標準出力への出力を results.txt にも同時に書き出すためのクラス"""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

if __name__ == "__main__":
    model = SAND_64_SAT_Model()
    # num_rounds:r, num_bounds:m
    r, b = 12, 58

    # 選択暗号文攻撃の立場: 12段目の差分を固定し, 暗号文側から平文側へ逆伝搬させる
    # 12-round differential distinguisher of SAND-64 (Table 3 of the paper)
    # ciphertext_diff = ("0492 9600", "0003 8000")  # (dX^12, dY^12)
    ciphertext_diff = None
    # 平文側の差分 (dX^0, dY^0) も固定する場合は次の行を有効にする
    # plaintext_diff = ("0046 0000", "00F0 0000")
    plaintext_diff = None

    # 実行結果をスクリプトと同じディレクトリの results.txt に出力する
    result_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.txt")

    with open(result_file, 'w') as f:
        stdout = sys.stdout
        sys.stdout = Tee(stdout, f)
        try:
            print(f"Running SAND-64 differential cryptanalysis (backward / chosen ciphertext):")
            print(f"Rounds: {r}, Bound: {b}")
            model.verify_ssb_model()
            print("-" * 50)
            model.run_without_key(r, b, ciphertext_diff, plaintext_diff)
            #model.run_backward_without_key(r, b, ciphertext_diff, plaintext_diff)
        finally:
            sys.stdout = stdout

    print("-" * 50)
    print(f"Results are saved in {result_file}")