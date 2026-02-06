# canfar_monitor.py

# Montior job submitted with the canfar client

import os
import sys

from requests.exceptions import ReadTimeout
from httpx import HTTPStatusError
from pydantic import ValidationError
from canfar.sessions import Session

import pandas as pd

from cs_util import args as cs_args
from cs_util import logging


class Log(object):

    def __init__(self):

        self.params_default()

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
            "kind": "headless",
            "status": "all",
            "destroy": False,
            "bulk": 1,
            "quiet": False,
        }

        self._short_options = {
            "kind": "-k",
            "status": "-s",
            "destroy": "-d",
            "bulk": "-b",
            "quiet": "-q",
        }

        self._types = {
            "destroy": "bool",
            "bulk": "int",
            "quiet": "bool",
        }

        self._help_strings = {
            "kind": "session kind, allowed are all, headless, notebook, ...; default is {}",
            "status": "procesing status, allowed are Pending, Running, Completed, all; default is {}",
            "destroy": "if True, destroy all displayed jobs; default is {}",
            "bulk": "destroy mode, allowed are 2, 1, 0; default is {}",
            "quiet": "quiet output",
        }

    def update_params(self):
        """Update Params.

        Update parameters.
        """
        pass

    def check_params(self):
        pass

    def print_info(self, df):

        if not self._params["quiet"]:
            # Prepare header and rows with left-justified formatting
            header = f"{'type':<10} {'status':<10} {'date':<12} {'time':<10} {'name':<20} {'id':<10}"
            print("=" * len(header))
            print(header)
            print("=" * len(header))

            for _, row in df.iterrows():
                print(f"{row['type']:<10} {row['status']:<10} {row['date']:<12} {row['time']:<10} {row['name']:<20} {row['id']:<10}")

            print("=" * len(header))

        print(f"{len(df)} jobs found")

    def filter(self):

        df = pd.DataFrame(self._info)

        if len(self._info) == 0:
            return df

        # Extract clean date/time
        try:
            # Remove 'T' and 'Z' to make it a proper datetime string
            dt_series = (
                df['startTime']
                    .str.replace('T', ' ', regex=False)
                    .str.replace('Z', '', regex=False)
                )
            # Extract date and time separately
            df['date'] = dt_series.str.split(' ').str[0]
            df['time'] = dt_series.str.split(' ').str[1]
        except KeyError:
            if self._params["verbose"]:
                print(
                    f"Column 'startTime' not found among {df.columns.tolist()},"
                    + f" setting dummy columns for date and time."
                )
            df['date'] = '1975-12-06'
            df['time'] = '12:04:00'
        except Exception as e:
            estr = f": {e}" if self._params["verbose"] else ""
            raise RuntimeError(
                f"An error occurred while processing dates: {estr}"
            ) from e

        # Select columns of interest
        cols = ['type', 'status', 'date', 'time', 'name', 'id']
        for col in cols:
            if not col in df:
                raise IndexError(f"Column {col} not in job info, which contains {df.columns.tolist()}") 
        df = df[['type', 'status', 'date', 'time', 'name', 'id']]

        # Get input status or None if missing
        status_filter = self._params.get("status", None)

        # If not None and not "all": filter
        if status_filter and status_filter != "all":
            df = df[df["status"] == status_filter]

        return df

    def run(self, args=None):

        if args is None:
            args = sys.argv

        obj = self
   
        obj.set_params_from_command_line(args)
        obj.update_params()
        obj.check_params()

        session = self.get_session()
        
        df_filtered = self.filter()

        self.print_info(df_filtered)

        if self._params["destroy"]:
            self.destroy(df_filtered, session)

    def get_kind(self):
        """Get Kind.

        Return kind for communication with canfar.
        In particular, returns None if kind in parametesr is "all".

        """
        return (
            self._params["kind"]
            if self._params["kind"] != "all"
            else None
        )

    def destroy(self, df, session):

        bulk = self._params["bulk"]

        if bulk == 2:
            kind = self.get_kind()

            if self._params["status"] != "all":
                status = self._params["status"]
            else:
                status = None

            try:
                print("Bulk destroy")
                result = session.destroy_with(
                    prefix="sp-",
                    kind=kind,
                    status=status,
                )
            except Exception as e:
                if self._params.get("verbose"):
                    print(f"Failed to destroy jobs: {e}")
        
            print("Success")
            return

        elif bulk == 1:
            # Recommended
            ids = df['id'].tolist()
            try:
                print("List destroy")
                result = session.destroy(ids)
                print(f"{len(ids)} sessions done")
            except Exception as e:
                estr = f": {e}" if self._params["verbose"] else ""
                print(f"failed{estr}")

        else:
            print("Loop destroy")
            for _, row in df.iterrows():
                session_id = row["id"]
                name = row.get("name", "<unknown>")
                print(f"Destroying {session_id} {name}", end=" ", flush=True)
                try:
                    result = session.destroy(session_id)
                    print(f"{result[session_id]} done")
                except Exception as e:
                    estr = f": {e}" if self._params["verbose"] else ""
                    print(f"failed{estr}")


    def get_session(self):

        # Initialize session manager
        session = Session()

        # Fetch information
        try:
            if not self._params["quiet"]:
                print(
                    f"Retrieving session info of kind={self._params['kind']},"
                    + f" status={self._params['status']}"
                )
            kind = self.get_kind()
            self._info = session.fetch(kind=kind)
            if len(self._info) == 0:
                print("No session fetched")
                return []
        except ReadTimeout as e:
            estr = f": {e}" if self._params["verbose"] else ""
            raise ReadTimeout(
                f"ReadTimeout occurred while fetching session info{estr}"
            ) from e
        except ValidationError as e:
            estr = f": {e}" if self._params["verbose"] else ""
            raise ValueError(
                f"Validation error occurred; incorrect kind{estr}"
            ) from e
        except HTTPStatusError as e:
            estr = f": {e}" if self._params["verbose"] else ""
            raise ValueError(
                f"HTTPStatusError occurred while fetching session info{estr}"
            ) from None
        except Exception as e:
            # Handle any other exceptions
            print(f"An unexpected error occurred: {e}")
            if self._params["verbose"]:
                raise
            else:
                sys.exit(2)
        else:
            # Runs only if no exception occurred
            if not self._params["quiet"]:
                print(f"Successfully fetched kind '{self._params['kind']}'")

        return session
