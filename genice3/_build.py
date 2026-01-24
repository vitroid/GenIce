"""
カスタムビルドフック: symlinkを実際のファイルにコピーする
"""
import os
import shutil
from pathlib import Path
from setuptools.command.build_py import build_py as _build_py


class BuildPy(_build_py):
    """symlinkを実際のファイルにコピーするカスタムビルドコマンド"""

    def copy_file(self, src, dst, preserve_mode=True, preserve_times=True, link=None, level=1):
        """ファイルをコピーする際に、symlinkの場合は実際のファイルをコピー"""
        # symlinkの場合は、ターゲットファイルをコピー
        if os.path.islink(src):
            # symlinkのターゲットを取得
            target = os.readlink(src)
            
            # 相対パスの場合は、srcのディレクトリからの相対パスとして解決
            if not os.path.isabs(target):
                src_dir = os.path.dirname(src)
                target_path = os.path.normpath(os.path.join(src_dir, target))
            else:
                target_path = target
            
            # ターゲットファイルが存在する場合のみコピー
            if os.path.exists(target_path):
                # ターゲットファイルをコピー（再帰的に呼び出し）
                return super().copy_file(
                    target_path, dst, preserve_mode, preserve_times, link, level
                )
            else:
                # ターゲットが見つからない場合は警告を出してスキップ
                print(f"Warning: symlink target not found: {src} -> {target_path}")
                return dst, False
        
        # 通常のファイルの場合は親クラスのメソッドを呼び出す
        return super().copy_file(src, dst, preserve_mode, preserve_times, link, level)




