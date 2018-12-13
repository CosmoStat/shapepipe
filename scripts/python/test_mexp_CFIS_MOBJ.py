# Test whether object ID's in object exposure files are unique.
# Expected output: (name, 0) for all files

from astropy.io import fits
import numpy as np
import glob
import os

path     = 'temp'
pattern  = 'CFIS_MOBJ'
obj_list = glob.glob('{}/{}*'.format(path, pattern))

for obj in obj_list:
    f = fits.open(obj)
    u, i = np.unique(f[1].data['ID'], return_counts=True)
    print(obj, len(f[1].data))
    for j in range(5):
        print(j, len(u[i == j]))



