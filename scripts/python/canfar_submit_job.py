#!/usr/bin/env python
# canfar_submit_job.py

# Submit job with the new canfar client (python)

import os
import sys
from canfar.sessions import Session
from datetime import datetime

from cs_util import args as cs_args
from cs_util import logging


class Job(object):

    def __init__(self):

        self.params_default()

        self._patch = os.environ["patch"] if "patch" in os.environ else "P0"


    def set_params_from_command_line(self, args):
        """Set Params From Command line.

        Only use when calling using python from command line.
        Does not work from ipython or jupyter.

        """
        # Read command line options
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def params_default(self):

        self._params = {
            "job": -1,
            "exclusive": None,
            "file_IDs": None,
            "psf": "psfex",
            "n": "0",
            "debug_out": None,
            "test": 0,
            "log": "log_job_tile_IDs.txt",
        }

        self._short_options = {
            "job": "-j",
            "exclusive": "-e",
            "file_IDs": "-f",
            "psf": "-p",
            "dry_run": "-n",
            "log": "-l",
        }

        self._types = {
            "test": "int",
        }

        self._help_strings = {
            "job": "Running JOB, bit-coded, default is {}",
            "exclusive": "Run exclusively given image",
            "file_IDs": "file containing IDs",
            "psf": "PSF model, allowed are 'psfex' and 'mccd'; default is {}",
            "dry_run": "dry run, no actual processing",
            "debug_out": "debug output file path, default:set automatically",
            "log": "output log file with job and tile IDs",
            "test": "test run level (2, 1, 0: no test, default)"
        }

    def update_params(self):
        """Update Params.

        Update parameters.
        """
        pass

    def check_params(self):

        if self._params["test"] == 2:
            return

        if self._params["job"] == "-1":
            raise ValueError("No job indicated, use option -j")
        if not self._params["exclusive"] and not self._params["file_ID"]:
            raise ValueError("No image ID(s) indicated, use option -e ID or -f file_IDs") 


    def set_options_base(self):

        opt = {}

        if self._params["test"] == 2:
            return ""

        if self._params["test"] == 1:
            opt["h"] = "--test"
        else:
            # Options
            opt["j"] = f"-j {self._params['job']}"
            opt["p"] = f"-p {self._params['psf']}"

            # Set identified for debug file name
            if self._params["test"] == 2:
                deb_ID = "test"
            elif self._params["exclusive"]:
                deb_ID = self._params["exclusive"]
            else:
                deb_ID = self._params["file_IDs"]

            cwd = os.getcwd()
            if self._params['debug_out']:
                opt["debut_out"] = f"--debug_out {self._params['debug_out']}"
            else:
                opt["debut_out"] = f"--debug_out {cwd}/debug/debug_{deb_ID}.txt"

            opt["d"] = f"-d {cwd}"
            opt["m_s"] = "-m 1 -s 1" if self._patch in ("P8", "P9") else ""

        options = ""
        for key in opt:
            options = f"{options} {opt[key]}"

        # Remove leading whitespace, problem in passing args below
        options = options.lstrip()

        return options

    def set_job_name(self, tile_ID):

        tile_ID_dash = tile_ID.replace(".", "-")

        job_name = f"sp-{self._patch}-{self._params['job']}-{tile_ID_dash}"

        if self._params["test"] == 1:
            job_name = f"{job_name}-test{test}"

        return job_name

    def set_command(self):

        if self._params["test"] == 2:
            cmd = f"{os.environ['HOME']}/test.sh"
        else:
            cmd = f"{os.environ['HOME']}/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"

        return cmd

    def run(self):

        # Initialize session manager
        session = Session()

        # Set job parameters
        version = "1.1"
        image = f"images.canfar.net/unions/shapepipe:{version}"

        options_base = self.set_options_base()

        cmd = self.set_command()

        if self._params["test"] == 2:
            tile_IDs = ["test"]
        elif self._params["exclusive"]:
            tile_IDs = [self._params["exclusive"]]
        else:
            with open(self._params["file_IDs"], "r") as f:
                tile_IDs = [line.strip() for line in f]

        print(f"Submitting {len(tile_IDs)} jobs to canfar queue")
        job_ids = []
        for tile_ID in tile_IDs:

            job_name = self.set_job_name(tile_ID)

            options = f"{options_base} -e {tile_ID}"

            print(f"Running '{cmd} {options}'")

            # Submit flexible job (default - auto-scaling)
            job_id = session.create(
                name=job_name,
                image=image,
                cmd=cmd,
                args=options,
            )

            print(f"Submitted job: {job_id}")
            job_ids.append(job_id)

        return job_ids, tile_IDs

    def write_job_tile_IDs(self, job_ids, tile_IDs):

        with open(self._params["log"], "w") as log:
            for job_id, tile_ID in zip(job_ids, tile_IDs):
                print(job_id, tile_ID, file=log)

def run_job(*args):

    # Create instance
    obj = Job()

    obj.set_params_from_command_line(args)
    obj.update_params()
    obj.check_params()

    job_ids, tile_IDs = obj.run()
    print(f"Submitting {len(job_ids)} jobs: done")

    obj.write_job_tile_IDs(job_ids, tile_IDs)

def main(argv=None):
    """Main

    Main program
    """
    if argv is None:
        argv = sys.argv[1:]
    run_job(*argv)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
