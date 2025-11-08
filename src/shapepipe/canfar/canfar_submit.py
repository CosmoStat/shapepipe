# canfar_submit_job.py

# Submit job with the canfar client

import os
import sys
import asyncio
from canfar.sessions import Session
from canfar.sessions import AsyncSession
from datetime import datetime

from cs_util import args as cs_args
from cs_util import logging


class Job(object):

    def __init__(self):

        self.params_default()

        self._patch = os.environ["patch"] if "patch" in os.environ else "P0"

        # Set job parameters
        version = "1.1"
        self._image = f"images.canfar.net/unions/shapepipe:{version}"


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
            "dry_run": 0,
            "test": 0,
            "log": "log_jobs.txt",
            "sync": "asnyc",
        }

        self._short_options = {
            "job": "-j",
            "exclusive": "-e",
            "file_IDs": "-f",
            "psf": "-p",
            "dry_run": "-n",
            "log": "-l",
            "sync": "-S",
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
            "test": "test run level (2, 1, 0: no test, default)",
            "sync": "job submission mode, allowed are async, sync; default is {}",
        }

    def update_params(self):
        """Update Params.

        Update parameters.
        """
        pass

    def check_params(self):

        if self._params["sync"] == "async" and not self._params["file_IDs"]:
            raise ValueError("asynchronous mode only possible when tile_IDs file given (-f)")

        if self._params["test"] >= 2:
            return

        if self._params["job"] == "-1":
            raise ValueError("No job indicated, use option -j")
        if not self._params["exclusive"] and not self._params["file_IDs"]:
            raise ValueError("No image ID(s) indicated, use option -e ID or -f file_IDs")

    def set_options_base(self):

        opt = {}

        cwd = os.getcwd()
        
        # asynchronous mode
        if self._params["sync"] == "async":
            opt["-f"] = f"-f {cwd}/{self._params['file_IDs']}"

        if self._params["test"] >= 2:
            opt["-dummy"] = "--dummy"
        else:
            if self._params["test"] == 1:
                opt["h"] = "--test"
            else:

                # Debug
                opt["h"] = "--test"

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
      
                if self._params['debug_out']:
                    opt["debut_out"] = f"--debug_out {cwd}/debug/{self._params['debug_out']}"
                else:
                    opt["debut_out"] = f"--debug_out {cwd}/debug/debug_{deb_ID}.txt"

                opt["d"] = f"-d {cwd}"
                opt["m_s"] = "-m 1 -s 1" if self._patch in ("P8", "P9") else ""
                if self._params["dry_run"] != 0:
                    opt["n"] = f"-n {self._params['dry_run']}"

        options = ""
        for key in opt:
            options = f"{options} {opt[key]}"

        self._options_base = options

    def set_job_name(self, suf=None):

        if suf is None:
            suf = self._params["exclusive"]

        suf1 = suf.replace(".", "-")
        
        job_name = f"sp-{self._patch}-{self._params['job']}-{suf1}"

        if self._params["test"] == 1:
            job_name = f"{job_name}-test-1"

        return job_name

    def set_command(self):

        if self._params["test"] >= 2:
            self._cmd = f"{os.environ['HOME']}/test/test.sh"
        elif self._params["sync"] == "async":
            self._cmd = f"{os.environ['HOME']}/shapepipe/scripts/sh/canfar_async_job.sh"
        else:
            self._cmd = f"{os.environ['HOME']}/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"

    def get_tile_IDs(self):

        if self._params["test"] >= 2:
            tile_IDs = [f"test-{i}" for i in range(self._params["test"])]
            suf = "test"
        elif self._params["exclusive"]:
            tile_IDs = [self._params["exclusive"]]
            suf = None
        else:
            with open(self._params["file_IDs"], "r") as f:
                tile_IDs = [line.strip() for line in f]
            suf = ""

        return tile_IDs, len(tile_IDs), suf

    async def run_async(self):
        async with AsyncSession() as session:

            _, n, suf = self.get_tile_IDs()

            print(f"Submitting {n} jobs to canfar queue")

            job_name = self.set_job_name(suf=suf)
            options = self._options_base.lstrip()
            print(f"Running async '[{self._cmd}] [{options}]'")
            
            try:
                sessions = await session.create(
                    name=job_name,
                    image=self._image,
                    cmd=self._cmd,
                    args=options,
                    replicas=n,
                )
                print("Success? sessions = ", sessions)
            except Exception as e:
                print(f"❌ CANFAR session.create() failed: {type(e).__name__}: {e}")
                raise

        return sessions

    def run_no_async(self):

        # Initialize session manager
        session = Session()

        tile_IDs, _, _ = self.get_tile_IDs()

        job_ids = []
        for tile_ID in tile_IDs:

            job_name = self.set_job_name(tile_ID)

            options = f"{self._options_base} -e {tile_ID}"

            # Remove leading whitespace, problem in passing args below
            options = options.lstrip()

            print(f"Running '[{self._cmd}] [{options}]'")

            # Submit flexible job (default - auto-scaling)
            job_id = session.create(
                name=job_name,
                image=self._image,
                cmd=self._cmd,
                args=options,
            )

            print(f"Submitted job: {job_id}")
            job_ids.append(job_id)

        return job_ids

    def write_job_tile_IDs(self, job_ids):

        tile_IDs, _, _ = self.get_tile_IDs()

        cwd = os.getcwd()
        log_path = f"{cwd}/{self._params['log']}"

        print(f"Writing job and tile ID log file {log_path}")
        with open(log_path, "w") as log:
            for job_id, tile_ID in zip(job_ids, tile_IDs):
                print(job_id, tile_ID, file=log)

    def run(self, args=None):

        obj = self

        if args is None:
            args = sys.argv
        obj.set_params_from_command_line(args)
        obj.update_params()
        obj.check_params()

        obj.set_options_base()
        obj.set_command()

        if obj._params["sync"] == "async":
            print("Async mode")
            job_ids = asyncio.run(obj.run_async())
        else:
            print("Sync mode")
            job_ids = obj.run_no_async()
        print(f"Submitting jobs: done")

        obj.write_job_tile_IDs(job_ids)
