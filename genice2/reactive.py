import inspect
import logging

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

logger = logging.getLogger(__name__)


class property_depending_on:
    """
    依存関係情報を保持し、遅延評価とキャッシュを行うディスクリプタ。
    """

    def __init__(self, *dependencies):
        # 依存関係の情報をインスタンス変数として保持する
        self._dependencies = dependencies
        self.func = None  # デコレートされる関数を保持

    def __call__(self, func):
        # デコレートされた関数を受け取り、それを保持し、自分自身（ディスクリプタインスタンス）を返す
        self.func = func
        return self

    def __set_name__(self, owner, name):
        # クラス定義時に、プロパティ名を記録
        self.attr_name = name

    def __get__(self, instance, owner=None):
        if instance is None:
            return self  # クラスからアクセスされた場合はディスクリプタ自身を返す

        # 2. キャッシュロジック（初回アクセス時のみ実行）
        if self.attr_name not in instance.__dict__:
            logger.debug(">>> %s を計算しています...", self.attr_name)
            # 関数を実行し、インスタンスの辞書に直接キャッシュする
            instance.__dict__[self.attr_name] = self.func(instance)

        return instance.__dict__[self.attr_name]


# ----------------------------------------------------
# 2. DependencyCacheMixin は変更なし (そのまま使用可能)
# ----------------------------------------------------


import inspect


class DependencyCacheMixin:
    """
    Observer Pattern を利用して、推移的キャッシュ無効化を実現する Mixin。
    """

    def __setattr__(self, name, value):
        # 1. 属性を更新する前に、以前の値を取得しておく（変更されたかどうかの判定用）
        old_value = getattr(self, name, object())

        # 2. 属性の更新を実行
        super().__setattr__(name, value)

        # 3. 値が実際に変更されたか、または依存元として定義されている属性かを確認
        #    （ここでは、簡略化のため、_で始まる属性が変更されたら通知をトリガーすると仮定）
        if name.startswith("_"):
            # numpy配列の場合は適切な比較を行う
            if HAS_NUMPY and isinstance(value, np.ndarray):
                if isinstance(old_value, np.ndarray):
                    value_changed = not np.array_equal(value, old_value)
                else:
                    # old_valueがnumpy配列でない場合は変更ありとみなす
                    value_changed = True
            elif HAS_NUMPY and isinstance(old_value, np.ndarray):
                # valueがnumpy配列でない場合は変更ありとみなす
                value_changed = True
            else:
                # 通常の比較（numpy配列でない場合）
                value_changed = value != old_value

            if value_changed:
                # 変更があった属性を「通知」の起点とする
                self._notify_dependents(name)

    def _notify_dependents(self, changed_attr_name):
        """
        指定された属性に依存する全てのプロパティのキャッシュを再帰的に無効化する。
        """
        # クラス内の全ての lazy_property (および標準 property) を走査
        for prop_name, prop in inspect.getmembers(
            self.__class__,
            lambda x: isinstance(x, (property, property_depending_on)),  # <<< 修正点
        ):

            # カスタム lazy_property であることと、依存関係があることを確認
            # lazy_property ディスクリプタは __init__ で _dependencies を持っている
            if hasattr(prop, "_dependencies"):
                # 依存関係のチェック:
                # 1. 直接マッチ（changed_attr_name が依存関係に含まれる）
                # 2. _で始まる属性の場合、_を削除した名前でもチェック
                #    （例: "_unitcell" -> "unitcell" もチェック）
                dependency_matched = changed_attr_name in prop._dependencies or (
                    changed_attr_name.startswith("_")
                    and changed_attr_name[1:] in prop._dependencies
                )

                if dependency_matched:
                    # 依存プロパティのキャッシュを無効化
                    if prop_name in self.__dict__:
                        logger.debug(
                            "--- 依存元 '%s' の変更により、キャッシュ '%s' を無効化 ---",
                            changed_attr_name,
                            prop_name,
                        )
                        del self.__dict__[prop_name]

                    # **重要: ここから連鎖が発生**
                    self._notify_dependents(prop_name)
                    # この再帰呼び出しが make のような連鎖的な無効化を実現します。


# lazy_property デコレータはそのまま使用
# @lazy_property('circumference') のように、直接的な依存関係を宣言する


def main():
    class Circle(DependencyCacheMixin):
        def __init__(self, radius):
            self._radius = radius  # 根源的な依存元

        @property
        def radius(self):
            return self._radius

        @radius.setter
        def radius(self, new_radius):
            if new_radius != self._radius:
                self._radius = new_radius  # __setattr__ がここでトリガーされる

        # 直接の依存先は _radius
        @property_depending_on("_radius")
        def circumference(self):
            print(">>> circumferenceを計算しています...")
            return 2 * 3.14159 * self.radius

        # 直接の依存先は circumference プロパティ
        @property_depending_on("circumference")
        def area(self):
            print(">>> areaを計算しています...")
            return self.circumference * self.radius / 2

    c = Circle(10)
    print(f"Area (1回目): {c.area}")
    # >>> circumferenceを計算しています...
    # >>> areaを計算しています...
    # Area (1回目): 314.159

    # 根源的な依存元を変更
    c.radius = 20

    # _radiusが変更され、以下の連鎖が発生:
    # 1. _radiusの変更通知 -> circumferenceのキャッシュ削除
    # 2. circumferenceのキャッシュ削除通知 -> areaのキャッシュ削除

    # 2回目アクセス時: 連鎖により両方再計算
    print(f"\nArea (2回目): {c.area}")
    # >>> circumferenceを計算しています...  <- 再計算
    # >>> areaを計算しています...          <- 再計算
    # Area (2回目): 1256.636


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
