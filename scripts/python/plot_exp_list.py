from astropy.io import ascii
from generic import plot
import pylab as plt

dat = ascii.read('exp_footprints.txt')
ax = plot.plot_init()
lw = 0.25


for d in dat:
    cx = [d['col2'], d['col4'], d['col6'], d['col8']]
    cy = [d['col3'], d['col5'], d['col7'], d['col9']]
    plt.plot(cx, cy, '-', linewidth=lw)
 

plt.savefig('exp_footprints.pdf')

