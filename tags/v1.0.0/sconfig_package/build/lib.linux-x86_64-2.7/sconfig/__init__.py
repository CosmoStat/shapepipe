"""!

@mainpage SConfig - Simple configuration file parser

@section features Main features:
- Key-Value pairs may belong to sections (as in MS Windows&trade; .ini files)
- Sub-sections
- Section inheritance
- Key values can be retrieved as: string, boolean, float, integer, list, tuple, dictionary, array, module
- Long line may be continued using the '\' character

See test/config.cfg file for examples.

@section Installation

- Download and unpack the file sconfig-1.0.0.tar.gz in the directory of your
  choice
- sconfig is a pure Python module and can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where sconfig-1.0.0.tar.gz was detared </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable </li>
</ol>

@section notes Release notes

- v1.0.0: initial stable version
"""

from sconfig import *
from sconfig_version import __version__
