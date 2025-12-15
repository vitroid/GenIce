from functools import wraps
import time
from logging import getLogger

# decorator


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
        logger.info(f"{func.__name__}: {msg.strip()}")
        try:
            return func(*args, **kwargs)
        finally:
            logger.info(f"{func.__name__}: end.")
    return _banner
