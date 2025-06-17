# -*- coding: utf-8 -*-
from importlib.metadata import version, PackageNotFoundError

def get_version():
    try:
        return version("genice2")
    except PackageNotFoundError:
        # インストール前の場合は、プロジェクトのルートディレクトリからバージョンを取得
        try:
            import toml

            with open("pyproject.toml") as f:
                t = toml.load(f)
                if "version" in t["project"]:
                    return t["project"]["version"]
                return t["tool"]["poetry"]["version"]
        except (ImportError, FileNotFoundError, KeyError):
            return "unknown"
