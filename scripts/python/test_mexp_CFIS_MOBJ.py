# Prints histograms of object ID's in object multi-exposure files.
# Expected output per file:
# (name, 0)
#  0 count_0
#  1 count_1
#  ...
# The largest counts should occur for n=3.

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



