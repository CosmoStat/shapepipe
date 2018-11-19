#!/usr/bin/env python

from astropy.io import fits

hdu_out = fits.HDUList()

for i in range(40):

    fname = 'cfisexp_flag-{:03d}-0.fits'.format(i)
    print(fname)

    hdu = fits.open(fname)
    hdu_out.append(hdu[0])

hdu_out.writeto('test.fits', overwrite=True)

