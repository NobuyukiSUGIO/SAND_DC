#!/usr/bin/env python3
"""Kissat の呼び出しと解の取り出し。ブロック長に依存しない。"""


class SolverMixin:
    """Solver 部分。変種固有のモデルクラスに mixin する。"""

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
