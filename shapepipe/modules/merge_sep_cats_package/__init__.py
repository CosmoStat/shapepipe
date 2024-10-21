r"""MERGE SEP CATS PACKAGE.

This package contains the module for ``merge_sep_cat``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: ``ngmix_runner``

:Input: SExtractor catalogue

:Output: SExtractor catalogue

Description
===========

This module merges separate input catalogues. This applies to
a run of the previous module that was split into several (parallel)
sub-processes to save wall-clock time. This is useful for time-consuming
modules, such as galaxy shape measurement e.g. via ``ngmix_runner``. This
module merges the outputs of those sub-processes back into a complete output
catalogue.

Note that the input catalogues should not contain identical objects.

Module-specific config file entries
===================================

N_SPLIT_MAX : int
    Number of input sub-catalogues
WARNING : str, optional
    Warning action if one or more input catalogue is missing; default is
    ``error`` (one of the actions defined
    `here <https://docs.python.org/3/library/warnings.html#warning-filter>`_)

"""

__all__ = ["merge_sep_cats.py"]
