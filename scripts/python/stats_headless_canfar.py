#!/usr/bin/env python

# Name: stats_headless_canfar.py

import sys
from skaha.session import Session


def main(argv=None):
    session = Session()

    n_headless = session.stats()["instances"]["headless"]

    print(n_headless)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
