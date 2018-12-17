#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii

from generic import plot

dat = ascii.read('exp_footprints.txt')
ax = plot.plot_init()
lw = 0.25

exp_names = np.unique(dat['name'])
colors = plot.color('tableau20')
col    = {}
for e in exp_names:
    col[e] = colors.next()

for d in dat:
    cx = [d['cx1'], d['cx2'], d['cx3'], d['cx4'], d['cx1']]
    cy = [d['cy1'], d['cy2'], d['cy3'], d['cy4'], d['cy1']]
    c  = col[d['name']]
    plt.plot(cx, cy, '-', color=c, linewidth=lw)
    plt.text(np.mean(cx), np.mean(cy), str(d['num']), color=c, fontsize=3,
             horizontalalignment='center', verticalalignment='center')
 

plt.savefig('exp_footprints.pdf')

