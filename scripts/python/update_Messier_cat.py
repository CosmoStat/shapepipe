#!/usr/bin/env python

from __future__ import print_function

import numpy as np
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord

# Download Messier catalogue
m_SB_cat    = Simbad.query_object("m [0-9]*", wildcard=True, verbose=True)

# Load exising, local Messier catalogue
cat_path = 'Messier_catalog.npy'
m_local_cat = np.load(cat_path)

# Loop over local catalogue
for m_local in m_local_cat:

    # Get Messier number
    no = m_local['No']

    # Get corresponding downloaded object
    m_SB  = m_SB_cat[no-1]

    # Downloaded coordinates
    sc_SB = SkyCoord(ra=str(m_SB['RA']), dec=str(m_SB['DEC']), unit=(u.hourangle, u.deg))

    print('{} {}     {} {}    '.format(no, m_SB['MAIN_ID'], m_local['ra'], m_local['dec']), end='')

    # Update local coordinates
    m_local['ra']  = sc_SB.ra.value
    m_local['dec'] = sc_SB.dec.value

    print('{} {}'.format(m_local['ra'], m_local['dec']))

cat_path = 'Messier_catalog_updated.npy'
np.save(cat_path, m_local_cat)

