"""!
@mainpage SCatalog - Simple catalog management

@section features Main features:
- Read and write text catalogs in tablulated text and SExtractor formats
- Read .FITS tables (the update table feature will be done in a later release) 

See test/sctest.py file for examples

@section Prerequisites 

- Python 2.6 or 2.7 on Unix or Mac OS
- AstroAsciiData v1.1 or greater (http://www.stecf.org/software/PYTHONtools/astroasciidata)
- pyfits v3.1.x

@section Installation

- Download and unpack the file scatalog-1.0.0.tar.gz in the directory of your
  choice
- SCatalog is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where scatalog-1.0.0.tar.gz was unpacked </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable </li>
  </ol>

@section notes Release notes

- v1.0.0: initial stable version
"""

from scatalog import *
from scatalog_version import __version__
