#!/usr/bin/env csh

# Script  remove_hdu0.sh
# Author: Martin Kilbinger
# Date:   02/2018

# Unzip .fits.fz files, write as .fits.
# Remove empty HDU #0 in CFIS weight images using cfitsio:imcopy

# Files .fits must not exist.


foreach i (*.r.weight.fits.fz)
	set new = `basename $i .fz`
	echo "imcopy $i\[1\] $new"
	imcopy $i\[1\] $new
end


