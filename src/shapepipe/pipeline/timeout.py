"""TIMEOUT.

This module defines a function for handling job timeout limits.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import signal
from functools import wraps


def with_timeout(timeout, log_file):
    """Timeout Limit Decorator.

    This method provides a timeout decorator for a given input function.

    Parameters
    ----------
    timeout : int
        Timeout limit in seconds
    log_file : logging.Logger
        Logging instance

    Raises
    ------
    TimeoutError
        For process exceeding timeout limit

    """
    def handler(signum, frame):
        raise TimeoutError(
            f'The process time exceeded {timeout}s in {log_file}'
        )

    def decorator(decorated):

        @wraps(decorated)
        def inner(*args, **kwargs):
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            return decorated(*args, **kwargs)

        return inner

    return decorator
