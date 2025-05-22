from functools import wraps
import time
from logging import getLogger
from typing import Any, Callable, Dict, TypeVar, cast
from icecream import ic

# decorator

T = TypeVar("T", bound=Callable[..., Any])


def debug_args(func: T) -> T:
    """関数の引数を表示するデコレータ"""

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        # 関数名を取得
        func_name = func.__name__

        ic(func_name)
        ic(args)
        ic(kwargs)

        # 元の関数を実行
        return func(*args, **kwargs)

    return cast(T, wrapper)


def timeit(func):
    @wraps(func)
    def _time_it(*args, **kwargs):
        logger = getLogger()
        start = int(round(time.time() * 1000))
        try:
            return func(*args, **kwargs)
        finally:
            end_ = int(round(time.time() * 1000)) - start
            logger.info(f"{func.__name__}: {end_ if end_ > 0 else 0} ms")

    return _time_it


# decorator


def banner(func):
    @wraps(func)
    def _banner(*args, **kwargs):
        logger = getLogger()
        lines = func.__doc__.splitlines()
        if len(lines[0]) == 0 and len(lines) > 1:
            msg = lines[1]
        else:
            msg = lines[0]
        # クラスメソッドの場合はクラス名のみを表示
        if hasattr(func, "__qualname__"):
            func_name = func.__qualname__.split(".")[0]
        else:
            func_name = func.__name__
        logger.info(f"{func_name}: {msg.strip()}")
        try:
            return func(*args, **kwargs)
        finally:
            logger.info(f"{func_name}: end.")

    return _banner
