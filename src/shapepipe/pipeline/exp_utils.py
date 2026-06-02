"""EXPOSURE UTILITIES.

Utility functions for accessing per-exposure runner outputs from the
tile level.  In the v2.0 pipeline, each single exposure is processed in
its own work directory ``<EXP_BASE_DIR>/<prefix>/<base>/``.  Tile-level
modules that need files produced by a per-exposure runner use the
functions here to discover those files by scanning the contributing
exposure directories.

:Author: Martin Kilbinger

"""

import glob
import os


def get_exp_output_files(
    exp_base_dir,
    exp_numbers_file,
    runner_name,
    file_pattern,
    file_ext,
    w_log=None,
):
    """Collect output files from a per-exposure runner for all tile exposures.

    For each exposure ID listed in ``exp_numbers_file``, the most recent
    output file matching ``<file_pattern>-<exp_base><file_ext>`` is located
    inside the exposure work directory::

        <exp_base_dir>/<exp_prefix>/<exp_base>/output/run_sp_*/<runner_name>/output/

    where ``exp_prefix = exp_id[:2]`` and ``exp_base`` is ``exp_id`` with the
    trailing letter stripped if present (e.g. ``2113864p`` → ``2113864``), or
    the full ``exp_id`` for numeric-only IDs (image simulations).

    Parameters
    ----------
    exp_base_dir : str
        Root directory that contains all per-exposure work directories,
        e.g. ``/arc/home/kilbinger/v2.0/exp``
    exp_numbers_file : str
        Path to a text file listing one exposure ID per line,
        e.g. ``exp_numbers-301-279.txt``
    runner_name : str
        Module runner directory name whose output is needed,
        e.g. ``split_exp_runner``
    file_pattern : str
        File name prefix used by the runner, e.g. ``headers``
    file_ext : str
        File extension including the leading dot, e.g. ``.npy``
    w_log : logging.Logger, optional
        Pipeline logger; ``None`` silences all logging

    Returns
    -------
    list of list
        ``[[path], [path], ...]`` — one single-element list per exposure,
        in the order they appear in ``exp_numbers_file``.  This format
        matches the ``input_file_list`` convention used throughout
        ShapePipe runners.

    Raises
    ------
    FileNotFoundError
        If ``exp_numbers_file`` does not exist, or if the expected output
        file is absent for one or more exposures.

    """
    if not os.path.exists(exp_numbers_file):
        raise FileNotFoundError(
            f"Exposure numbers file not found: {exp_numbers_file}"
        )

    with open(exp_numbers_file) as fh:
        exp_ids = [line.strip() for line in fh if line.strip()]

    if w_log:
        w_log.info(
            f"Collecting {runner_name}/{file_pattern}*{file_ext} "
            f"for {len(exp_ids)} exposures from {exp_base_dir}"
        )

    file_list = []
    missing = []

    for exp_id in exp_ids:
        # Directory structure mirrors run_job_canfar_v2.0.sh:
        #   exp_prefix = first 2 chars of exp_id  (e.g. "21")
        #   exp_base   = exp_id without trailing letter if present (e.g. "2113864"),
        #                or full exp_id for numeric-only ids (image sims)
        exp_prefix = exp_id[:2]
        exp_base = exp_id[:-1] if exp_id[-1].isalpha() else exp_id

        pattern = os.path.join(
            exp_base_dir,
            exp_prefix,
            exp_base,
            "output",
            "run_sp_*",
            runner_name,
            "output",
            f"{file_pattern}-{exp_base}{file_ext}",
        )
        matches = sorted(glob.glob(pattern))

        if matches:
            chosen = matches[-1]  # most recent run (lexicographic sort on datetime)
            file_list.append([chosen])
            if w_log:
                w_log.debug(f"  {exp_id}: {chosen}")
        else:
            missing.append(exp_id)
            if w_log:
                w_log.warning(f"  {exp_id}: no match for {pattern}")

    if missing:
        raise FileNotFoundError(
            f"No {runner_name} output found for "
            f"{len(missing)} exposure(s): {missing}"
        )

    if w_log:
        w_log.info(f"Found {len(file_list)} exposure output files")

    return file_list


def get_exp_output_dirs(
    exp_base_dir,
    exp_numbers_file,
    runner_name,
    w_log=None,
):
    """Return the most-recent output directory of runner_name for each exposure.

    For each exposure ID listed in ``exp_numbers_file``, the most recent
    runner output directory is located at::

        <exp_base_dir>/<exp_prefix>/<exp_base>/output/run_sp_*/<runner_name>/output/

    Parameters
    ----------
    exp_base_dir : str
        Root directory containing all per-exposure work directories.
    exp_numbers_file : str
        Path to a text file listing one exposure ID per line.
    runner_name : str
        Module runner directory name whose output directory is needed,
        e.g. ``psfex_runner``.
    w_log : logging.Logger, optional
        Pipeline logger.

    Returns
    -------
    list of str
        One output directory path per exposure, in the order they appear
        in ``exp_numbers_file``.

    Raises
    ------
    FileNotFoundError
        If ``exp_numbers_file`` does not exist or a runner output directory
        is absent for one or more exposures.

    """
    if not os.path.exists(exp_numbers_file):
        raise FileNotFoundError(
            f"Exposure numbers file not found: {exp_numbers_file}"
        )

    with open(exp_numbers_file) as fh:
        exp_ids = [line.strip() for line in fh if line.strip()]

    if w_log:
        w_log.info(
            f"Collecting {runner_name}/output dirs "
            f"for {len(exp_ids)} exposures from {exp_base_dir}"
        )

    dir_list = []
    missing = []

    for exp_id in exp_ids:
        exp_prefix = exp_id[:2]
        exp_base = exp_id[:-1] if exp_id[-1].isalpha() else exp_id

        pattern = os.path.join(
            exp_base_dir,
            exp_prefix,
            exp_base,
            "output",
            "run_sp_*",
            runner_name,
            "output",
        )
        matches = sorted(glob.glob(pattern))

        if matches:
            chosen = matches[-1]  # most recent run (lexicographic sort on datetime)
            dir_list.append(chosen)
            if w_log:
                w_log.debug(f"  {exp_id}: {chosen}")
        else:
            missing.append(exp_id)
            if w_log:
                w_log.warning(f"  {exp_id}: no match for {pattern}")

    if missing:
        raise FileNotFoundError(
            f"No {runner_name} output directory found for "
            f"{len(missing)} exposure(s): {missing}"
        )

    if w_log:
        w_log.info(f"Found {len(dir_list)} exposure output directories")

    return dir_list
