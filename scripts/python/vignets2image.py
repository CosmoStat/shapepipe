#!/usr/bin/env python

from shapepipe.pipeline import file_io as io

from astropy.io import fits

import re
import os
import sys
import glob
import numpy as np


def map_vignet(img_arr, dtype):
    """Map vignet
    Map vignet on one single image.

    Parameters
    ----------
    img_arr : numpy.ndarray
        Array of vignets to map
    dtype : str
        dtype of the data

    Returns
    -------
    img_map : numpy.ndarray
        Array containing all the vignets mapped on one single image
    nx : int
        Number of objects along one side (assumed square image)
    """

    n_obj = img_arr.shape[0]
    xs = img_arr[0].shape[0]
    ys = img_arr[0].shape[1]

    nx = int(np.sqrt(n_obj))
    if nx*nx != n_obj:
        nx += 1
    ny = nx

    img_map=np.ones((xs*nx,ys*ny), dtype=dtype)

    ii=0
    jj=0
    for i in range(n_obj):
        if jj>nx-1:
            jj=0
            ii+=1
        img_map[ii*xs:(ii+1)*xs,jj*ys:(jj+1)*ys]=img_arr[i]
        jj+=1

    return img_map, nx



def main(argv=None):

    path = '.'
    pattern = 'UNIONS_'
    prefix_out = 'img_'

    files = glob.glob('{}/{}*.fits'.format(path, pattern))

    for input_path in files:
        output_path = '{}{}'.format(prefix_out, os.path.basename(input_path))

        hdu = fits.open(input_path)

        image, nx = map_vignet(hdu[2].data['VIGNET'], 'float32')
        print('file = {}, nx = {}'.format(input_path, nx))

        fout = io.FITSCatalog(output_path, open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        fout.save_as_fits(image, image=True)

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

