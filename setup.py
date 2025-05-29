from setuptools import setup
import os
import sys

# カレントディレクトリをPythonパスに追加
sys.path.insert(0, os.path.abspath("."))

from genice2.build import BuildPyCommand

setup(
    cmdclass={"build_py": BuildPyCommand},
)
