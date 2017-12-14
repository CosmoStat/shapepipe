"""!
@mainpage SCatalog - Simple catalog management

@section features Main features:
- Read and write text catalogs in tablulated text and SExtractor formats
- Read .FITS tables (the update table feature will be done in a later release) 

See test/sctest.py file for examples

@section Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- AstroAsciiData v1.1 or greater (http://www.stecf.org/software/PYTHONtools/astroasciidata)
- Astromatics astromatic_wrapper v1.0 (github.com/fred3m/astromatic_wrapper)
- Astropy v1.3.2 or greater (http://www.astropy.org/)

@section Installation

- Download and unpack the file scatalog-2.0.1.tar.gz in the directory of your
  choice
- SCatalog is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where scatalog-2.0.1.tar.gz was unpacked </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable </li>
  </ol>

@section notes Release notes

- v1.0.0: initial stable version
- v1.0.2: reimplementation of FITS interface using fitsio instead of pyfits
- v2.0.0: reimplementation of FITS catalogs with the astropy.io.fits package;
          added LDACFITSCatalog class for minimal support for Astromatics FITS LDAC format
- v2.0.1: fix some bug;
          FITSCatalog can handle "every" fits format including SExtractor_LDAC fits format 
"""

from scatalog import *
from scatalog_version import __version__
