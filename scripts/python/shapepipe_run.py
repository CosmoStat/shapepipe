#!/usr/bin/env python

"""SHAPEPIPE RUN SCRIPT

This script runs the shape measurement pipeline.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from sys import exit

from shapepipe.run import run


def main(args=None):

    run(args)


if __name__ == "__main__":
    exit(main())
