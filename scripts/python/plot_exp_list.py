#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

from generic import plot

dat = ascii.read('exp_footprints.txt')
ax = plot.plot_init()
lw = 0.25


col = 'g'
for d in dat:
    cx = [d['col2'], d['col4'], d['col6'], d['col8'], d['col2']]
    cy = [d['col3'], d['col5'], d['col7'], d['col9'], d['col3']]
    plt.plot(cx, cy, '{}-'.format(col), linewidth=lw)
    plt.text(np.mean(cx), np.mean(cy), str(d['col1']), color=col, fontsize=3,
             horizontalalignment='center', verticalalignment='center')
 

plt.savefig('exp_footprints.pdf')

