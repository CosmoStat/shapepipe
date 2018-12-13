# Test whether object ID's in object exposure files are unique.
# Expected output: (name, 0) for all files

from astropy.io import fits
import numpy as np
import glob
import os

path_hdu = '{}/data/hdu'.format(os.environ['HOME'])
pattern  = 'cfisexp-obj'
obj_list = glob.glob('{}/{}*'.format(path_hdu, pattern))

for obj in obj_list:
    f = fits.open(obj)
    u, i = np.unique(f[2].data['ID'], return_inverse=True)
    print(obj, len(u[np.bincount(i) > 1]))



