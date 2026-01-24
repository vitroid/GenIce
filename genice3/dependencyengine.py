import inspect
from logging import getLogger
import time


class DependencyEngine:
    logger = getLogger("DependencyEngine")

    def __init__(self):
        self.registry = {}  # { 'output_name': function }
        self.cache = {}

    def task(self, func):
        """デコレータ: 関数名 = 生成される変数名 として登録"""
        self.registry[func.__name__] = func
        return func

    def resolve(self, target: str, inputs: dict):
        """ゴール(target)から逆算して計算する"""

        # 1. 既に計算済み or 入力として与えられているならそれを返す
        if target in inputs:
            return inputs[target]
        if target in self.cache:
            return self.cache[target]

        # 2. 生成ルール（関数）を探す
        if target not in self.registry:
            raise ValueError(f"Don't know how to make '{target}'")

        func = self.registry[target]

        # 3. その関数の引数（依存先）を調べて、再帰的に解決する
        sig = inspect.signature(func)
        dependencies = {}
        for param_name in sig.parameters:
            dependencies[param_name] = self.resolve(param_name, inputs)

        # 4. 実行して結果を保存
        self.logger.info(f"Executing: {target}")  # 実行ログ
        now = time.time()
        result = func(**dependencies)
        delta = time.time() - now
        self.logger.info(f"{delta:.4f} sec for {target}")
        self.cache[target] = result
        return result


def main():
    # --- ライブラリ定義 ---
    engine = DependencyEngine()

    @engine.task
    def reaction_A(raw_material: int):
        return raw_material * 2

    @engine.task
    def reaction_B(reaction_A: int, catalyst: int):
        # reaction_A が必要だと自動判別される
        return reaction_A + catalyst

    # --- ユーザー利用 ---
    # ユーザーは「reaction_B が欲しい」と言うだけ
    result = engine.resolve(
        target="reaction_B", inputs={"raw_material": 10, "catalyst": 5}
    )
    print(result)


if __name__ == "__main__":
    main()






