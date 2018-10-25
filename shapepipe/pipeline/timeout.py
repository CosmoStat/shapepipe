# -*- coding: utf-8 -*-

"""TIMEOUT

This module defines a function for handling job timeout limits.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import multiprocessing as mp
from functools import wraps


def with_timeout(timeout, log_file):
    """ Timeout Limit Decorator

    This method provides a timeout decorator for a given input function.

    Parameters
    ----------
    timeout : int
        Timeout limit in seconds
    log_file : Logger
        Logging instance

    Raises
    ------
    mp.TimeoutError
        For process exceeding timeout limit

    Notes
    -----
    Thise method was taken from: https://github.com/joblib/joblib/pull/366

    """

    def decorator(decorated):

        @wraps(decorated)
        def inner(*args, **kwargs):
            pool = mp.pool.ThreadPool(1)
            async_result = pool.apply_async(decorated, args, kwargs)
            try:
                return async_result.get(timeout)
            except mp.TimeoutError:
                raise TimeoutError('The process time exceeded {}s in '
                                   '{}'.format(timeout, log_file))
        return inner

    return decorator
