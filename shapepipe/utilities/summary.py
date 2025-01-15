"""SUMMARY

Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import sys
import os
import re
import fnmatch

import logging

from collections import Counter

from tqdm import tqdm

print("summaary v1.1")


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

        #print(inds_special)
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
    path_main: str, optional
        main (left-most) part of output directory, default is "."       
    path_left: str, optional
        left (first) part of output directory, default is "./output"
    output_subdirs: str, optional
        output subdirectories if not `None`; default is `None`
    path_right: str, optional
        right (last) part of output subdir suffix if not `None`;
        default is `None`
    output_path_missing_IDs: list, optional
        output path of missing ID, if `None` (default) will be
        given by job bit and module.
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
        path_main=".",
        path_left="output",
        output_subdirs=None,
        path_right=None,
        output_path_missing_IDs=None,
        verbose=False,
    ):
        self._bit = bit
        self._run_dir = set_as_list(item=run_dir, n=len(modules))
        self._modules = modules
        self._key_expected = set_as_list(item=key_expected, n=len(modules))
        self._n_mult = set_as_list(item=n_mult, n=len(modules))
        self._pattern = set_as_list(item=pattern, n=len(modules), default="")
        self._path_main = path_main
        self._path_left = path_left
        self._output_subdirs = output_subdirs or [""]
        self._path_right = set_as_list(
            path_right, len(modules), default="."
        )
        self._output_path_missing_IDs=output_path_missing_IDs
        self._verbose = verbose

    def print_intro(self):
        """Print Intro.

        Print header line for job statistics.

        """
        logging.info(f" # Job {self._bit}:")

    @classmethod
    def print_stats_header(self):
        """Print Stats Header.

        Print overall header information for stats output.

        """
        logging.info(
            "module                          expected     found   miss_expl"
            + " missing uniq_miss  fr_found"
        )
        logging.info("=" * 100)

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
    def replace_dot_dash(self, numbers):

        results = [re.sub("\.", "-", number) for number in numbers]

        return results

    @classmethod
    def replace_dash_dot_if_tile(self, numbers):

        pattern = re.compile(r"(\d{3})-(\d{3})")
        results = [pattern.sub(r"\1.\2", number) for number in numbers]
    
        return results

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

    @classmethod
    def write_IDs_to_file(self, output_path, IDs):
        """Write IDs to file.

        Write list if image IDs to text file.

        Parameters
        ----------
        output_path: str
            output file path
        IDs: list
            image IDs

        """
        IDs_dot = self.replace_dash_dot_if_tile(IDs)
        with open(output_path, "w") as f_out:
            for ID in IDs_dot:
                print(ID, file=f_out)

    def output_missing(
        self,
        module,
        idx,
        par_runtime=None,
    ):
        """Output Missing.

        Writes IDs of missing images to disk.

        """
        key_expected = self._key_expected[idx]
        names_in_dir = self._names_in_dir[idx]
        paths_in_dir = self._paths_in_dir[idx]
        n_mult = self._n_mult[idx]

        list_expected = get_par_runtime(par_runtime, key_expected, kind="list")
        
        # Count image IDs in names that were found earlier

        ## Extract image IDs from names
        IDs = []
        pattern = re.compile(
            r"(?:\d{3}-\d{3}|\d{7}-\d+|\d{7})"
        )
        for name, path in zip(names_in_dir, paths_in_dir):
            match = pattern.search(name)
            if match:
                IDs.append(match.group())
            else:
                raise ValueError(f"No ID found in {name}")


        ## Count occurences
        ID_counts = Counter(IDs)

        ## Add to missing if ocurence less than n_mult
        missing_IDs = []
        for ID in list_expected:
            if ID_counts[ID] < n_mult:
                missing_IDs.append(ID)

        n_all = len(missing_IDs)
        missing_IDs_unique = self.get_unique(missing_IDs)
        n_unique = len(missing_IDs_unique)

        if n_unique > 0:
            if not self._output_path_missing_IDs:
                output_path = (
                    f"{self._path_main}/summary/missing_job_{self._bit}"
                    + f"_{module}.txt"
                )
            else:
                output_path = self._output_path_missing_IDs[idx]
            #print("MKDEBUG", missing_IDs_unique)
            self.write_IDs_to_file(output_path, missing_IDs_unique)

        return missing_IDs_unique

    def output_missing_job(self):
        output_path = (
            f"{self._path_main}/summary/missing_job_{self._bit}_all.txt"
        )

        missing_IDs_all = set(self._missing_IDs_job)

        if len(missing_IDs_all) > 0:
            self.write_IDs_to_file(output_path, missing_IDs_all)
        else:
            #logging.warning("no missing IDs in output_missing_job")
            if os.path.exists(output_path):
                os.unlink(output_path)

    @classmethod
    def get_last_full_path(self, base_and_subdir, matches):
        """Get Last Full Path

        Return full path of last file in list.

        """
        # Sort according to creation time
        matches_sorted = sorted(
            matches,
            key=lambda entry: entry.name,
        )

        # Get most recent one
        last = matches_sorted[-1]

        # Get full path
        full_path = os.path.join(base_and_subdir, last.name)

        return full_path

    @classmethod
    def get_module_output_dir(self, full_path, module):
        """Get Module Output Dir.

        Return output directory name for given module.

        """
        directory = f"{full_path}/{module}/output"

        # Some modules have special requirements 
        if module == "setools_runner":
            directory = f"{directory}/rand_split"

        return directory

    def get_matches_final(self, directory, idx):

        # Loop over files
        # os.path.whether exists is twice faster than try/except

        if os.path.exists(directory):
            pattern =  f"{self._pattern[idx]}*"
            for entry2 in os.scandir(directory):
                if (
                    entry2.is_file()
                    and (
                        fnmatch.fnmatch(entry2.name, pattern)
                    )
                    and entry2.stat().st_size > 0
                ):
                    # Append matching files
                    self._names_in_dir[idx].append(entry2.name)
                    self._paths_in_dir[idx].append(
                        os.path.join(directory, entry2.name)
                    )

    def get_names_in_dir(self, iterable, module, idx):

        # Initialise output file names and paths
        self._names_in_dir[idx] = []
        self._paths_in_dir[idx] = []

        # Loop over subdirs
        for jdx, subdir in enumerate(iterable):
            base_and_subdir = (
                f"{self._path_main}/"
                + f"{self._path_left}/{subdir}/"
                + f"{self._path_right[idx]}"
            )
            if self._verbose:
                print(f"**** base_and_subdir {base_and_subdir}")

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

                    if self._verbose:
                        print("**** Matching entries: ", end="")
                        for match in matches:
                            print(match.name)

                    full_path = self.get_last_full_path(
                        base_and_subdir, matches
                    )

                    # Get module output directory
                    directory = self.get_module_output_dir(
                        full_path, module
                    )
                    if self._verbose:
                        print(f"**** Output dir = {directory}")

                    # Find matching file names and paths
                    self.get_matches_final(directory, idx)
            else:
                if self._verbose:
                    print(f"Directory {base_and_subdir} not found")

    def update_subdirs(self, par_runtime):
        """Update Subdirs.

        Update subdir names with runtime information if required.

        """
        if not isinstance(self._output_subdirs, list):
            self._output_subdirs = get_par_runtime(
                par_runtime, self._output_subdirs, kind="list"
            )

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
        self.update_subdirs(par_runtime)

        # Initialise variables
        self._names_in_dir = {}
        self._paths_in_dir = {}
        self._missing_IDs_job = []
        n_missing_job = 0

        # Loop over modules
        for idx, module in enumerate(self._modules):
            if indices is not None and idx not in indices:
                continue

            if self._verbose:
                print(f"** module {module}")

            # Look over subdirs
            iterable = self._output_subdirs
            if len(iterable) > 1 and self._verbose:
                iterable = tqdm(iterable, desc="subdirs", leave=False)

            if self._verbose:
                print(f"*** subdirs {self._output_subdirs}")

            # Get output file names and paths
            self.get_names_in_dir(
                iterable,
                module,
                idx,
            )

            # If expected is string: Update parameter with runtime value
            # and set as integer
            if isinstance(self._key_expected[idx], str):
                n_expected_base = get_par_runtime(
                    par_runtime, self._key_expected[idx], kind="n"
                )
            else:
                n_expected_base = self._key_expected[idx]

            # Get some numbers
            n_found = len(self._names_in_dir[idx])
            n_expected = n_expected_base * self._n_mult[idx]
            n_missing = n_expected - n_found

            n_missing_explained = 0
            if n_missing > 0:
                # TODO: make check_special class function, deal with
                # paths, names in dir
                if False and module == "setools_runner":
                    (
                        self._paths_in_dir[idx],
                        self._names_in_dir[idx],
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
                missing_IDs = self.output_missing(
                    module,
                    idx,
                    par_runtime=par_runtime,
                )
                n_missing_job += n_missing
                self._missing_IDs_job.extend(missing_IDs)

        # Empty line after job
        logging.info("")

        # Write missing IDs for entire job to file
        #if n_missing_job > 0:
        self.output_missing_job()


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
    if (
        combined_key == "list_3*n_shdus+n_exposures"
        and combined_key not in par_runtime
    ):
        print("{combined_key} not set, TBD")
        return []

    return par_runtime[combined_key]


def print_par_runtime(par_runtime, verbose=True):
    # Print runtime parameter values
    if True:
        logging.info("")
        logging.info("===========")
        logging.info("par_runtime")
        logging.info("-----------")
        for key, value in par_runtime.items():
            if not key.startswith("list"):
                logging.info(f"{key:30s} {value:6d}")
            else:
                logging.info(f"{key:29s} {len(value):6d} entries")
        logging.info("===========")
        logging.info("")
