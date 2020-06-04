from astropy.io import fits
import numpy as np
import os
import re
from tqdm import tqdm

l = os.listdir('.')

d_tmp = fits.getdata(l[0], 1)                                                                
#new_dt = np.dtype(d_tmp.dtype.descr + [('TILE_ID', '>i4')])
#d = np.zeros(d_tmp.shape, dtype=new_dt)
d = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
for key in d_tmp.dtype.names:
    d[key] = d_tmp[key]
d['TILE_ID'].fill(int(''.join(re.findall('\d+', l[0]))))
for i in tqdm(l[1:], total=len(l)-1):
    if ('final_cat' not in i) | ('.npy' in i):
        continue
    d_tmp = fits.getdata(i, 1)
    #new_dt = np.dtype(d_tmp.dtype.descr + [('TILE_ID', '>i4')])
    #dd = np.zeros(d_tmp.shape, dtype=new_dt)
    dd = np.zeros(d_tmp.shape, dtype=d_tmp.dtype)
    for key in d_tmp.dtype.names:
        dd[key] = d_tmp[key]
    dd['TILE_ID'].fill(int(''.join(re.findall('\d+', i))))
    d = np.concatenate((d, dd))

np.save('final_cat.npy', d)