# -*- coding: utf-8 -*-

"""CREATE LOG EXP HEADER 

This script merge the "headers" file output of the split_exp_runner.py module.
It create a binnary file that contain the wcs of each CCDs for each exposures. 

:Author: Axel Guinot

"""

import glob
import os
import sys
import re

import numpy as np


def main(input_path, output_name):
    """
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