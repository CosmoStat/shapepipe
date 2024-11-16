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

print("summaary v1.4")


def init_par_runtime(list_tile_IDs):

    # Numbers updated at runtime
    par_runtime = {}

    par_runtime["n_tile_IDs"] = len(list_tile_IDs)
    par_runtime["list_tile_IDs"] = list_tile_IDs

    return par_runtime


def update_par_runtime_after_find_exp(par_runtime, all_exposures):

    # Single-exposure images
    par_runtime["n_exposures"] = len(all_exposures)
    par_runtime["list_exposures"] = all_exposures

    # Single-HDU single exposure images
    n_CCD = 40
    par_runtime["n_shdus"] = get_par_runtime(par_runtime, "exposures") * n_CCD
    par_runtime["list_shdus"] = get_all_shdus(all_exposures, n_CCD)

    return par_runtime


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


def check_special_one(module, path):

    with open(path) as f_in:
        lines = f_in.readlines()
        for line in lines:
            entry = line.rstrip()

            if module == "setools_runner":
                m = re.search("Nb stars = (\S*)", line)
                if m:
                    value = int(m[1])
                    if value < 2:
                        code = 0
                        msg = (
                            f"Not enough stars for random split:"
                            + f"  #stars = {value}"
                        )
                        return msg, code
                    break
                m = re.search("Mode computation failed", line)
                if m:
                    code = 1
                    msg = "Mode computation of stellar locus failed"
                    return msg, code

            if module == "psfex_interp_runner":
                m = re.search("Key N_EPOCH not found", line)
                if m:
                    code = 2
                    msg = "N_EPOCH not in SEx cat, rerun job 16"
                    return msg, code

                m = re.search(
                    "ValueError: cannot reshape array of size 0 into shape",
                    line,
                )
                if m:
                    code = 3
                    msg = "found array of size 0"
                    return msg, code


    return None, None 


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
    path_output: str, optional
        module output path, default is "output"
    output_path_missing_IDs: list, optional
        output path of missing ID, if `None` (default) will be
        given by job bit and module.
    special: bool, optional
        if True check output file content for special messages;
        default is False
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
        path_output="output",
        output_path_missing_IDs=None,
        special=False,
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
            path_right,
            len(modules),
            default=".",
        )
        self._path_output = set_as_list(
            path_output,
            len(modules),
            default="output",
        )
        self._output_path_missing_IDs=output_path_missing_IDs
        self._special = set_as_list(
            special,
            len(modules),
            default=False,
        )
        self._path_right = set_as_list(path_right, len(modules), default=".")
        self._output_path_missing_IDs = output_path_missing_IDs
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
            "module                          expected     found"
            + "   missing uniq_miss  fr_found"
        )
        logging.info("=" * 100)

    def print_stats(
        self,
        module,
        n_expected,
        n_found,
        n_special,
        n_missing,
        idx,
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
        n_special: int
            number of special cases
        n_missing: int
            number of missing files
        idx: int
            module index

        """
        module_str = module
        
        if not self._special[idx]:
            if n_expected > 0:
                fraction_found = n_found / n_expected
            else:
                fraction_found = 1
                
            n_missing_per_mult = n_missing / self._n_mult[idx]
            
        else:
            module_str = f"{module_str} (special)"
            n_found = n_special
            n_missing = -1
            n_missing_per_mult = -1
            fraction_found = n_found / n_expected
            n_expected = -1

        logging.info(
            f"{module_str:30s} {n_expected:9d} {n_found:9d}"
            + f" {n_missing:9d}"
            + f" {n_missing_per_mult:9.1f} {fraction_found:9.1%}"
        )

    @classmethod
    def is_ID_in_str(self, ID, path):
        if ID in path:
            return True

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
        if len(IDs_dot) > 0:
            # Write IDs to file
            with open(output_path, "w") as f_out:
                for ID in IDs_dot:
                    print(ID, file=f_out)
        elif os.path.exists(output_path):
            # Remove preivous obsolete ID file
            os.unlink(output_path)

    def check_special(self, module, idx):

        messages = {}
        
        if self._special[idx]:
            
            # Loop over input file names and paths
            for name, path in zip(self._names_in_dir[idx], self._paths_in_dir[idx]):

                # Check if special case is found
                msg, code = check_special_one(module, path)
                if msg:
                    # First time occurance: create empty list for this code 
                    if code not in messages:
                        messages[code] = [msg]
                    else:
                        # Append file name, message, and code 
                        messages[code].append(f"{name} {code} {msg}")

            if len(messages) > 0:
                # Loop over codes = key in messages dict
                for code in messages:
                    # Create output file for this code
                    output_path = (
                        f"{self._path_main}/summary/special_job_{self._bit}"
                        + f"_{module}_{code}.txt"
                    )
                    # Write all messages
                    with open(output_path, "w") as f_out:
                        for msg in messages[code]:
                            print(msg, file=f_out)

        # Count all special cases = sum of cases over all codes 
        n_all = sum([len(messages[code]) for code in messages])
        return n_all
        
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

        # Get file name pattern
        if module != "split_exp_runner" or self._bit != "2":
            pattern = re.compile(r"(?:\d{3}-\d{3}|\d{7}-\d+|\d{7})")
        else:
            # split_exp_runner with sp_local=0: input is exp, output is shdu
            # (images) and exp (header); ignore hdu number.
            # If sp_local=1 set bit to != 2
            pattern = re.compile(
                r"(?:\d{3}-\d{3}|\d{7})"
            )

        ## Extract image IDs from names
        IDs = []
        for name, path in zip(names_in_dir, paths_in_dir):

            match = pattern.search(name)
            if match:
                ID = match.group()
                IDs.append(ID)
            else:
                msg = f"No ID found in {name}"
                #raise ValueError(msg)
                print(f"Warning: msg, continuing")

        # For split_exp_runner P8, IDs now contain exps and sdus,
        # not matching mult.

        ## Count occurences
        ID_counts = Counter(IDs)

        ## Add to missing if ocurence less than n_mult
        missing_IDs = []
        for ID in list_expected:
            if ID_counts[ID] < n_mult:
                missing_IDs.append(ID)

        n_all = len(missing_IDs)
        missing_IDs_unique = self.get_unique(missing_IDs)

        if not self._output_path_missing_IDs:
            # Default name using bit and module
            output_path = (
                f"{self._path_main}/summary/missing_job_{self._bit}"
                + f"_{module}.txt"
            )
        else:
            # User-defined name (e.g. ngmix_runner_X)
            output_path = self._output_path_missing_IDs[idx]
        self.write_IDs_to_file(output_path, missing_IDs_unique)

        return missing_IDs_unique

    def output_missing_job(self):
        output_path = (
            f"{self._path_main}/summary/missing_job_{self._bit}_all.txt"
        )

        missing_IDs_all = set(self._missing_IDs_job)

        self.write_IDs_to_file(output_path, missing_IDs_all)

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
    def get_module_output_dir(self, full_path, module, path_output):
        """Get Module Output Dir.

        Return output directory name for given module.

        """
        directory = f"{full_path}/{module}/{path_output}"

        return directory

    def get_matches_final(self, directory, idx):

        # Loop over files
        # os.path.whether exists is twice faster than try/except

        if os.path.exists(directory):
            pattern = f"{self._pattern[idx]}*"
            for entry2 in os.scandir(directory):
                if (
                    entry2.is_file()
                    and (fnmatch.fnmatch(entry2.name, pattern))
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
                        base_and_subdir,
                        matches,
                    )

                    # Get module output directory
                    directory = self.get_module_output_dir(
                        full_path,
                        module,
                        self._path_output[idx],
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

            n_special = self.check_special(
                module,
                idx,
            )

            # Print statistics
            self.print_stats(
                module,
                n_expected,
                n_found,
                n_special,
                n_missing,
                idx,
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
        # if n_missing_job > 0:
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
                # logging.info(f"{key:30s} {len(value):6d} entries")
                pass
        logging.info("===========")
        logging.info("")
