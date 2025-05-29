from setuptools.command.build_py import build_py
import os
import shutil
from logging import getLogger


def clone_to(dst=".genice2/genice2", delete=False):
    print(f"Cloning to {dst}")
    if delete:
        print(f"Deleting {dst}")
        shutil.rmtree(dst, ignore_errors=True)

    # 親ディレクトリを作成
    os.makedirs(os.path.dirname(dst), exist_ok=True)

    def copy_file(src, dst):
        print(f"Copying {src} to {dst}")
        if os.path.islink(src):
            src = os.path.realpath(src)
            print(f"  Following symlink to {src}")
        shutil.copy(src, dst)

    def copy_dir(src_dir, dst_dir):
        print(f"Copying directory {src_dir} to {dst_dir}")
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            dst_path = os.path.join(dst_dir, item)
            if os.path.isdir(src_path):
                copy_dir(src_path, dst_path)
            else:
                copy_file(src_path, dst_path)

    # genice2ディレクトリをコピー
    copy_dir("genice2", dst)

    # デバッグ用：コピー後のディレクトリ内容を表示
    print(f"Contents of {dst}:")
    print(os.listdir(dst))


class BuildPyCommand(build_py):
    def run(self):
        print("Running BuildPyCommand.run()")
        # .genice2ディレクトリを作成し、genice2/の中の全ファイルを.genice2/genice2/にコピー
        clone_to(delete=True)

        # dstの中身を表示し、/Users/matto/に保存。
        with open("/Users/matto/dst.txt", "w") as f:
            f.write(os.listdir(dst))
        # ビルド先を指定
        self.build_lib = ".genice2"

        # パッケージのディレクトリ構造を設定
        self.package_dir = {"": ".genice2"}
        self.packages = [
            "genice2",
            "genice2.cli",
            "genice2.formats",
            "genice2.lattices",
            "genice2.groups",
            "genice2.loaders",
            "genice2.molecules",
        ]

        # 通常のビルド処理を実行
        super().run()


if __name__ == "__main__":
    clone_to(delete=True)
