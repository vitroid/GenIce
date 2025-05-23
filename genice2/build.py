from setuptools.command.build_py import build_py
import os
import shutil


def clone_to(dst=".genice2"):
    if not os.path.exists(dst):
        os.makedirs(dst)

    def copy_file(src, dst):
        if os.path.islink(src):
            src = os.path.realpath(src)
        shutil.copy(src, dst)

    def copy_dir(src_dir, dst_dir):
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            dst_path = os.path.join(dst_dir, item)
            if os.path.isdir(src_path):
                copy_dir(src_path, dst_path)
            else:
                copy_file(src_path, dst_path)

    copy_dir("genice2", dst)


class BuildPyCommand(build_py):
    def run(self):
        # make prepareを実行
        # subprocess.run(["make", "prepare"], check=True)

        # .genice2ディレクトリを作成し、genice2/の中の全ファイルを.genice2/genice2/にコピー。ただし、シンボリックリンクはリンク先のファイルに置きかえる。

        # 通常のビルド処理を実行
        super().run()


if __name__ == "__main__":
    clone_to()
