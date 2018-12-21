#!/usr/bin/env python

import os

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

pl = {}
for d in dat:
    cx = [d['cx1'], d['cx2'], d['cx3'], d['cx4'], d['cx1']]
    cy = [d['cy1'], d['cy2'], d['cy3'], d['cy4'], d['cy1']]
    c  = col[d['name']]

    # Plot CCD box
    pl[d['name']], = plt.plot(cx, cy, '-', color=c, linewidth=lw)

    # Write CCD/HDU number
    plt.text(np.mean(cx), np.mean(cy), str(d['num']), color=c, fontsize=3,
             horizontalalignment='center', verticalalignment='center')

# Legend
labels = []
plots  = []
for e in exp_names:
    labels.append(os.path.basename(e))
    plots.append(pl[e])
legend = plt.legend(plots, labels, loc=0, fontsize='x-small')
plt.gca().add_artist(legend)


plt.xlabel('R.A. [deg]')
plt.ylabel('DEC [deg]')

plt.savefig('exp_footprints.pdf')

