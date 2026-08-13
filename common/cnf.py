#!/usr/bin/env python3
"""CNF 構築の基本部品と CNF ファイルの出力。ブロック長に依存しない。"""


class CnfMixin:
    """Cnf 部分。変種固有のモデルクラスに mixin する。"""

    def rotate_left(self, data, bits):
        """Rotate a list to the left by a given number of bits."""
        return data[bits:] + data[:bits]

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
