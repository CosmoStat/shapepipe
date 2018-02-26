#!/usr/local/bin/csh

# Script  remove_hdu0.sh
# Author: Martin Kilbinger
# Date:   02/2018

# Remove empty HDU #0 in CFIS weight images using cfitsio:imcopy


foreach i (*.r.weight.fits.fz)
	set new = `basename $i .fz`
	echo "~/share/cexamples/imcopy $i\[1\] $new"
	#~/share/cexamples/imcopy $i\[1\] $new
end


