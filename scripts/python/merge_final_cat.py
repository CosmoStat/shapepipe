#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script merge_final_cat.py

Merge all final catalogues, created by ShapePipe module 'make_catalogue_runner',
into a joined numpy binary file.

:Authors: Axel Guinot

:Date: 2020

:Package: ShapePipe
"""

from astropy.io import fits
import numpy as np
import os
import sys
import re
from tqdm import tqdm


def main(argv=None):

    try:
        path = argv[1]
    except:
        path = '.'

    l = os.listdir(path=path)
    lpath = []
    for this_l in l:
        lpath.append(os.path.join(path, this_l))

    # Determine number of columns and keys
    d_tmp = fits.getdata(lpath[0], 1)
    d = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
    for key in d_tmp.dtype.names:
        d[key] = d_tmp[key]

    #new_dt = np.dtype(d_tmp.dtype.descr + [('TILE_ID', '>i4')])
    #d = np.zeros(d_tmp.shape, dtype=new_dt)

    d['TILE_ID'].fill(int(''.join(re.findall('\d+', l[0]))))

    # Read all final catalogues and merge
    for i in tqdm(lpath[1:], total=len(lpath)-1):
        if ('final_cat' not in i) | ('.npy' in i):
            continue

        #new_dt = np.dtype(d_tmp.dtype.descr + [('TILE_ID', '>i4')])
        #dd = np.zeros(d_tmp.shape, dtype=new_dt)

        d_tmp = fits.getdata(i, 1)
        dd = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
        for key in d_tmp.dtype.names:
            dd[key] = d_tmp[key]
        dd['TILE_ID'].fill(int(''.join(re.findall('\d+', i))))
        d = np.concatenate((d, dd))

    # Save merged catalogue as numpy binary file
    np.save('final_cat.npy', d)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

