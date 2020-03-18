# -*- coding: utf-8 -*-

"""CREATE LOG EXP HEADER

This script merges the "headers" file output of the split_exp_runner.py module.
It creates a binnary file that contain the WCS of each CCD for each exposure.

:Author: Axel Guinot

"""

import glob
import os
import sys
import re

import numpy as np


def main(input_path, output_name):
    """ Main

    Merge the individual files into a main one.

    Parameters
    ----------
    input_path : str
        Path to the input directory.
    output_name : str
        Path to the output file.
        
    """

    l = glob.glob(input_path + '/header*')

    d = {re.split('headers-', os.path.splitext(os.path.split(ll)[1])[0])[1]: np.load(ll) for ll in l}

    np.save(output_name, d)


if __name__ == '__main__':

    argv = sys.argv

    try:
        input_path = argv[1]
    except:
        input_path = '.'

    try:
        output_path = argv[2]
    except:
        output_path = '.'

    output_name = output_path + '/log_exp_headers.npy'

    main(input_path, output_name)
