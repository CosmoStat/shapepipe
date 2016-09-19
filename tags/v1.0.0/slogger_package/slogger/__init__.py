"""!
@mainpage FileLogger - Simple and lightweight file logging facility

@section Installation

- Download and unpack the file slogger-1.0.0.tar.gz in the directory of your
  choice
- This is a pure Python module and it can then be installed using the standard
  procedure: 

  <ol>
   <li> create a console or terminal session </li>
   <li> cd to the directory where slogger-1.0.0.tar.gz was detared </li>
   <li> enter the command: @code sudo python setup.py install @endcode </li>
   <li> if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable </li>
</ol> 
"""

from file_logger import *
from version import __version__

