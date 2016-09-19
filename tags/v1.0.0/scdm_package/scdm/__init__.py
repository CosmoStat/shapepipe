"""!
@mainpage SCDM - Simple Coordinate Descent Minimizer

Author: Marc Gentile (marc.gentile@epfl.ch)

@section Description

SCDM implements a variant of the Coordinate Descent Algorithm 
(http://en.wikipedia.org/wiki/Coordinate_descent).

@section Prerequisites 

- Mandatory:
   - Python 2.6 or 2.7 on Unix or Mac OS
   - slogger v1.x: logging module

@section Installation

Download and unpack the file scdm-0.6.0.tar.gz in the directory of your
choice. The SCDM module is a pure Python module and can thus be installed using the standard
procedure: 

<ol>
<li> create a console or terminal session </li>
<li> cd to the directory where scdm-0.6.0.tar.gz was unpacked </li>
<li> enter the command: command: @code sudo python setup.py install @endcode </li>

<li> if no root access to th default installation direcory, 
      use: @code sudo python setup install --home=target_directory 
           @endcode to install the module to target_directory. This directory must be
           included in the @c PYTHONPATH environment variable </li>
</ol> 

"""

from scdm import *
from scdm_version import __version__
