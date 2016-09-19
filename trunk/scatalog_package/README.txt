SCatalog v1.0.2
===============

SCatalog is a simple catalog management module.

Main features:
- Read and write text catalogs in tablulated text and SExtractor formats
- Read .FITS tables (the update table feature will be done in a later release)

See test/sctest.cfg file for examples

Prerequisites 
-------------

- Python 2.6 or 2.7 on Unix or Mac OS
- AstroAsciiData v1.1 or greater (http://www.stecf.org/software/PYTHONtools/astroasciidata)
- fitsio 0.9.7 or greater (https://github.com/esheldon/fitsio)

Installation
------------

- Download and unpack the file scatalog-1.0.2.tar.gz in the directory of your
  choice
- SCatalog is a pure Python module and can then be installed using the standard
  procedure: 

  - create a console or terminal session 
  - cd to the directory where scatalog-1.0.2.tar.gz was unpacked 
  - enter the command: @code sudo python setup.py install @endcode 
  - if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable 
 
Release notes
-------------

- v1.0.0: initial stable version
- v1.0.2: reimplementation of FITS interface using fitsio instead of pyfits
