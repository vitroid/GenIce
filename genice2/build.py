import subprocess
from setuptools.command.build_py import build_py


class BuildPyCommand(build_py):
    def run(self):
        # make prepareを実行
        subprocess.run(["make", "prepare"], check=True)
        # 通常のビルド処理を実行
        super().run()
