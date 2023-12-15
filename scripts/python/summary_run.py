#!/usr/bin/env python

import sys
import os
import re

import logging

from collections import Counter

from tqdm import tqdm


def get_IDs_from_file(path):
    """Get IDs From File.

    Return IDs from text file. Removes letters and replaces
    dots "." with dashes "-".

    Parameters
    ----------
    path: str
        input file path

    Returns
    --------
    list
        IDs

    """
    numbers = []
    with open(path) as f_in:
        for line in f_in:
            entry = line.rstrip()
            number = re.sub("[a-zA-Z]", "", entry)
            numbers.append(number)

    return numbers


def replace_dot_dash(numbers):

    results = [re.sub("\.", "-", number) for number in numbers]

    return results


def replace_dash_dot_if_tile(numbers):

    pattern = re.compile(r"(?<!\d{3}-\d{3})-")
    results = [pattern.sub(".", number) for number in numbers]
    
    return results


def replace_dot_dash(numbers):

    results = [re.sub("\.", "-", number) for number in numbers]

    return results



def get_all_exposures(exp_number_file_list, verbose=False):
    """Get All Exposures.

    Return all exposure names from a list of text files.

    Parameters
    ----------
    exp_number_list: list
        input file names

    """
    exposures = set()
    for idx, path in enumerate(exp_number_file_list):
        exps = get_IDs_from_file(path)
        exposures.update(exps)

    return list(exposures)


def get_all_shdus(exposures, n_CCD):
    """Get All SHDUs.

    Return all single-exposure single-HDU (CCD) IDs.

    Parameters
    ----------
    exposures: list
        exposure names
    n_CCD: int
        number of CCDs per exposure

    Returns
    --------
    list
        single-exposure single-HDU IDs

    """
    shdus = []
    for exposure in exposures:
        for idx_CCD in range(n_CCD):
            shdus.append(f"{exposure}-{idx_CCD}")

    return shdus


def set_as_list(item=None, n=None, default=1):
    """Set As List.

    Return input as list.

    Parameters
    -----------
    item: str, int, or list, optional
        input item(s); default is None, in which
        case the return is [1] * n
    n: int, optional
        number of list items to return, default is None,
        in which case the number will be set to 1. If item and
        n are not None, n has to be equal to len(item)
    default: int, optional
        value to return if item is not given;
        default is 1

    Raises
    -------
    IndexError
        if n != len(item)

    Returns
    -------
    list
        input item(s) as list
    """
    my_n = n or 1

    if not item:
        result = [default] * my_n
    elif not isinstance(item, list):
        result = [item] * my_n
    else:
        result = item
        if len(item) != my_n:
            raise IndexError(f"item has length {len(item)} != {n}")

    return result


def check_special(module, paths_in_dir, names_in_dir):

    if module == "setools_runner":
        inds_special = []
        for idx in range(len(paths_in_dir)):
            base_path = paths_in_dir[idx].replace(names_in_dir[idx], "")

            stats_dir = f"{base_path}/../stat"
            stats_files = os.listdir(stats_dir)
            if len(stats_files) != 1:
                raise ValueError(
                    f"Expected exactly one stats file in {stats_dir}, not"
                    + f" {len(stats_files)}"
                )

            stats_path = os.path.join(stats_dir, stats_files[0])
            with open(stats_path) as f_in:
                lines = f_in.readlines()
                for line in lines:
                    entry = line.rstrip()
                    m = re.search(line, "Nb stars = (\S*)")
                    if m:
                        value = int(m[2])
                        if value == 0:
                            inds_special.append(idx)
                        else:
                            print(f"b stars = {value}, not special")
                        break

        print(inds_special)
        for idx in inds_special:
            paths_in_dir.pop(idx)
            names_in_dir.pop(idx)

        return paths_in_dir, names_in_dir, len(inds_special)


class job_data(object):
    """Job Data.

    Class to handle a job.

    Parameters
    ----------
    bit: int
        bit-coded job number
    run_dir: str or list
        run directory(ies)
    modules: list
        module names
    key_expected: int or str
        number of expected output files; if str: will be updated
        with runtime value
    n_mult: int or list, optional
        multiplicity of output files, default `None`, in which
        case it is set to 1
    pattern: list, optional
        if not None, file pattern to match; defafult is `None`
    output_dir: str, optional
        output directory, defaul is "./output"
    output_subdirs: str, optional
        output subdirectories if not `None`; default is `None`
    output_subdirs_suffix: str, optional
        output subdir suffix if not `None`; default is `None`
    verbose: bool, optional
        verbose output if True; default is False

    """
    def __init__(
        self,
        bit,
        run_dir,
        modules,
        key_expected,
        n_mult=None,
        pattern=None,
        output_dir="./output",
        output_subdirs=None,
        output_subdirs_suffix=None,
        verbose=False,
    ):
        self._bit = bit
        self._run_dir = set_as_list(item=run_dir, n=len(modules))
        self._modules = modules
        self._key_expected = set_as_list(item=key_expected, n=len(modules))
        self._n_mult = set_as_list(item=n_mult, n=len(modules))
        self._pattern = set_as_list(item=pattern, n=len(modules), default="")
        self._output_dir = output_dir
        self._output_subdirs = output_subdirs or [""]
        self._output_subdirs_suffix = set_as_list(
            output_subdirs_suffix, len(modules), default="."
        )
        self._verbose = verbose

    def print_intro(self):
        """Print Intro.

        Print header line for job statistics.

        """
        logging.info(f" (Job {self._bit})")

    @classmethod
    def print_stats_header(self):
        """Print Stats Header.

        Print overall header information for stats output.

        """
        logging.info(
            "module                          expected     found   miss_expl"
            + " missing uniq_miss  fr_found"
        )
        logging.info("=" * 50)

    @classmethod
    def print_stats(
        self,
        module,
        n_expected,
        n_found,
        n_missing_explained,
        n_missing,
        n_mult,
    ):
        """Print Stats.

        Print output file statistics.

        Parameters
        ----------
        module: str
            module name
        n_expected: int
            number of expected files
        n_found: int
            number of found files
        n_missing_explained: int
            number of missing but explained files
        n_missing: int
            number of missing files
        n_mult: int
            multipicity

        """
        if n_expected > 0:
            fraction_found = n_found / n_expected
        else:
            fraction_found = 1

        n_missing_per_mult = n_missing / n_mult

        logging.info(
            f"{module:30s} {n_expected:9d} {n_found:9d}"
            + f" {n_missing_explained:9d} {n_missing:9d}"
            + f" {n_missing_per_mult:9.1f} {fraction_found:9.1%}"
        )

    @classmethod
    def is_ID_in_str(self, ID, path):
        if ID in path:
            return True
        #if re.sub("\.", "-", ID) in path:
            #return True
        #return False

    @classmethod
    def is_not_in_any(self, ID, list_str):
        return not any(ID in string for string in list_str)

    @classmethod
    def get_unique(self, names):
        n_all = len(names)
        names_unique = list(set(names))
        n_unique = len(names_unique)

        if n_all != n_unique:
            if True:  # self._verbose:
                logging.warning(
                    f"{n_all - n_unique} duplicates removed from {n_all} IDs"
                )

        return names_unique

    def output_missing(
        self, module, key_expected, names_in_dir, n_mult, par_runtime=None
    ):
        output_path = f"missing_job_{self._bit}_{module}.txt"

        list_expected = get_par_runtime(par_runtime, key_expected, kind="list")

        pattern = re.compile(r"\d+[\d-]+")
        IDs = []
        for name in names_in_dir:
            match = pattern.search(name)
            if match:
                IDs.append(match.group())
            else:
                raise ValueError(f"No ID found in {name}")

        ID_counts = Counter(IDs)
        missing_IDs = [
            ID for ID, count in ID_counts.items()
            if count < n_mult
        ]

        n_all = len(missing_IDs)
        missing_IDs_unique = self.get_unique(missing_IDs)
        n_unique = len(missing_IDs_unique)

        if n_unique > 0:
            IDs_dot = replace_dash_dot_if_tile(missing_IDs_unique)
            with open(output_path, "w") as f_out:
                for ID in IDs_dot:
                    print(ID, file=f_out)

        return missing_IDs_unique

    def output_missing_job(self, missing_IDs):
        output_path = f"missing_job_{self._bit}_all.txt"

        missing_IDs_all = set(missing_IDs)

        if len(missing_IDs_all) > 0:
            with open(output_path, "w") as f_out:
                for ID in missing_IDs_all:
                    print(ID, file=f_out)
        else:
            logging.warning("no missing IDs in output_missing_job")

    def get_names_in_dir(self, iterable, module, idx):

        # Initialise output file names
        names_in_dir = []
        paths_in_dir = []

        # Loop over subdirs
        for jdx, subdir in enumerate(iterable):
            base_and_subdir = (
                f"{self._output_dir}/{subdir}/"
                + f"{self._output_subdirs_suffix[idx]}"
            )
            if os.path.isdir(base_and_subdir):

                matches = []

                # Loop over entries (files and dirs)
                with os.scandir(base_and_subdir) as entries:
                    for entry in entries:

                        # Append directory name if matches module
                        if entry.name.startswith(self._run_dir[idx]):
                            matches.append(entry)

                    # This entry does not match module -> next
                    if not matches:
                        continue

                    # Sort according to creation time
                    matches_sorted = sorted(
                        matches,
                        key=lambda entry: entry.name,
                    )

                    # Get most recent one
                    last = matches_sorted[-1]

                    # Get full path
                    full_path = os.path.join(base_and_subdir, last.name)

                    # Get module output directory
                    directory = f"{full_path}/{module}/output"

                    # Some modules have special requirements 
                    if module == "setools_runner":
                        directory = f"{directory}/rand_split"

                    #if os.path.exists(directory):
                    try:
                        with os.scandir(directory) as entries2:
                            # if entry2.is_file()
                            files = [
                                entry2.name
                                for entry2 in entries2
                                if entry2.name.startswith(self._pattern[idx])
                            ]
                            names_in_dir.extend(files)
                            paths_in_dir.extend(
                                [os.path.join(directory, file)
                                for file in files]
                            )
                    except FileNotFoundError:
                        pass
                    except Exeption as e:
                        print(f"Unknown error {e}")

        return names_in_dir, paths_in_dir

    def check_numbers(self, par_runtime=None, indices=None):
        """Check Numbers.

        Check output file numbers and IDs.

        Parameters
        ----------
        par_runtime : dict, optional
            runtime parameter. default is None
        indices: list, optional
            if not None (default), only check modules corresponding
            to indices

        """
        # Update subdirs if not already set as list
        if not isinstance(self._output_subdirs, list):
            self._output_subdirs = get_par_runtime(
                par_runtime, self._output_subdirs, kind="list"
            )

        self._paths_in_dir = {}
        self._missing_IDs_unique = []

        n_missing_job = 0

        # Loop over modules
        for idx, module in enumerate(self._modules):
            if indices is not None and idx not in indices:
                continue

            # Look over subdirs
            iterable = self._output_subdirs
            if len(iterable) > 1 and self._verbose:
                iterable = tqdm(iterable, desc="subdirs", leave=True)

            # Get output file names and paths
            names_in_dir, paths_in_dir = self.get_names_in_dir(
                iterable,
                module,
                idx,
            )

            self._paths_in_dir[idx] = paths_in_dir

            # If expected is string: Update parameter with runtime value
            # and set as integer
            if isinstance(self._key_expected[idx], str):
                n_expected_base = get_par_runtime(
                    par_runtime, self._key_expected[idx], kind="n"
                )
            else:
                n_expected_base = self._key_expected[idx]

            # Get some numbers
            n_found = len(names_in_dir)
            n_expected = n_expected_base * self._n_mult[idx]
            n_missing = n_expected - n_found

            n_missing_explained = 0
            if False and n_missing > 0:
                if module == "setools_runner":
                    (
                        paths_in_dir,
                        names_in_dir,
                        n_missing_explained,
                    ) = check_special(module, paths_in_dir, names_in_dir)

                n_missing = n_missing - n_missing_explained

            # Print statistics
            self.print_stats(
                module,
                n_expected,
                n_found,
                n_missing_explained,
                n_missing,
                self._n_mult[idx],
            )

            # Write missing IDs for module to file
            if n_missing > 0:
                missing_IDs_unique = self.output_missing(
                    module,
                    self._key_expected[idx],
                    names_in_dir,
                    self._n_mult[idx],
                    par_runtime=par_runtime,
                )
                n_missing_job += n_missing
                self._missing_IDs_unique.extend(missing_IDs_unique)

        # Write missing IDs for entire job to file
        if n_missing_job > 0:
            self.output_missing_job(self._missing_IDs_unique)


def get_par_runtime(par_runtime, key, kind="n"):
    """Get Par RunTime.

    Return runtime parameter value.

    Parameters
    ----------
    par_runtime: dict
        runtime parameter
    key: str
        key

    """
    combined_key = f"{kind}_{key}"
    return par_runtime[combined_key]


def print_par_runtime(par_runtime, verbose=True):
    # Print runtime parameter values
    if verbose:
        logging.info("")
        logging.info("===========")
        logging.info("par_runtime")
        logging.info("-----------")
        for key, value in par_runtime.items():
            if not key.startswith("list"):
                logging.info(f"{key:30s} {value:6d}")
            else:
                logging.info(f"{key:29s} [{len(value):6d}]")
        logging.info("===========")
        logging.info("")


def main(argv=None):
    # Set default parameters
    # p_def = params_default()

    # Command line options
    # options, args = parse_options(p_def)

    # if check_options(options) is False:
    # return 1

    # param = update_param(p_def, options)

    verbose = True
    log_file_name = "summary_log.txt"
    handlers = [logging.FileHandler(log_file_name), logging.StreamHandler()]
    logging.basicConfig(
        level=logging.INFO, format="%(message)s", handlers=handlers
    )

    main_dir = "."
    retrieve = "vos"
    tile_ID_path = f"{main_dir}/tile_numbers.txt"

    # tile IDs with dots
    list_tile_IDs_dot = get_IDs_from_file(tile_ID_path)

    # tile IDs with dashes
    list_tile_IDs = replace_dot_dash(list_tile_IDs_dot)
    n_tile_IDs = len(list_tile_IDs)
    n_CCD = 40

    par_runtime = {}

    par_runtime["n_tile_IDs"] = n_tile_IDs
    par_runtime["list_tile_IDs"] = list_tile_IDs

    jobs = {}

    if retrieve == "vos":
        n_link = 2
    else:
        n_link = 1

    jobs["1"] = job_data(
        1,
        "run_sp_GitFeGie_",
        [
            "get_images_runner_run_1",
            "find_exposures_runner",
            "get_images_runner_run_2",
        ],
        ["tile_IDs", "tile_IDs", "exposures"],
        pattern=["CFIS_", "", ""],
        n_mult=[1 * n_link, 1, 3],
        output_dir=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["2"] = job_data(
        2,
        ["run_sp_Uz", "run_sp_exp_SpMh", "run_sp_exp_SpMh_2023-12"],
        ["uncompress_fits_runner", "merge_headers_runner", "split_exp_runner"],
        ["tile_IDs", 0, "3*n_shdus+n_exposures"],
        n_mult=[1, 1, 1],
        output_dir=f"{main_dir}/output",
        verbose=verbose,
    )

    # TODO: rename run dirs to run_sp_tile_Ma, run_sp_exp_Ma
    jobs["4"] = job_data(
        4,
        "run_sp_combined_flag",
        ["mask_runner_run_1"],
        ["tile_IDs"],
        output_dir=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["8"] = job_data(
        8,
        "run_sp_combined_flag",
        ["mask_runner_run_2"],
        ["shdus"],
        output_dir=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["16"] = job_data(
        16,
        "run_sp_tile_Sx",
        ["sextractor_runner"],
        ["tile_IDs"],
        n_mult=2,
        output_dir=f"{main_dir}/tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    # TODO setools_runner output/mask
    jobs["32"] = job_data(
        32,
        [
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
        ],  # "run_sp_exp_Pi"],
        [
            "sextractor_runner",
            "setools_runner",
            "psfex_runner",
        ],  # "psfex_interp_runner"],
        "shdus",
        n_mult=[2, 2, 2],  # 1],
        output_dir=f"{main_dir}/exp_runs",
        output_subdirs="shdus",
        output_subdirs_suffix="output",
        verbose=verbose,
    )

    jobs["64"] = job_data(
        "64",
        "run_sp_tile_PsViSmVi",
        [
            "psfex_interp_runner",
            "vignetmaker_runner_run_1",
            "spread_model_runner",
            "vignetmaker_runner_run_2",
        ],
        "tile_IDs",
        n_mult=[1, 1, 1, 4],
        output_dir=f"{main_dir}/tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    job_data.print_stats_header()

    for key in "1":
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime, indices=[0, 1])

        all_exposures = get_all_exposures(job._paths_in_dir[1], verbose=True)
        par_runtime["n_exposures"] = len(all_exposures)

        job.check_numbers(par_runtime, indices=[2])

    # Update runtime parameter
    par_runtime["n_shdus"] = get_par_runtime(par_runtime, "exposures") * n_CCD
    par_runtime["n_3*n_shdus+n_exposures"] = 3 * get_par_runtime(
        par_runtime, "shdus"
    ) + get_par_runtime(par_runtime, "exposures")
    par_runtime["list_shdus"] = get_all_shdus(all_exposures, n_CCD)

    print_par_runtime(par_runtime, verbose=verbose)

    #for key in ["2", "4", "8", "16", "32", "64"]:
    for key in ["32", "64"]:
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
