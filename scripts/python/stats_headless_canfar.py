#!/usr/bin/env python

# Name: stats_headless_canfar.py

# Caution: Does not show all running or pending
# headless jobs, for some reason.

import sys
from skaha.session import Session


def main(argv=None):

    print(
        "# Depreciated, does not show pending jobs; use stats_jobs_canfar.sh",
        file=sys.stderr,
    )

    session = Session()

    n_headless = session.stats()["instances"]["headless"]

    print(n_headless)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
