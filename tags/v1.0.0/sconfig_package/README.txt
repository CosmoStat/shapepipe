SConfig v1.0.0
==============

SConfig is a simple configuration file parser. 

Main features:
- Key-Value pairs may belong to sections (as in MS Windows' .ini files)
- Sub-sections
- Section inheritance
- Key values can be retrieved as: string, boolean, float, integer, list, dictionary
- long line may be continued using the '\' character

See test/config.cfg file for examples.

Installation

- Download and unpack the file sconfig-1.0.0.tar.gz in the directory of your
  choice
- sconfig is a pure Python module and can then be installed using the standard
  procedure: 

  - create a console or terminal session 
  - cd to the directory where sconfig-1.0.0.tar.gz was detared 
  - enter the command: @code python sudo setup.py install @endcode 
  - if no root access, use: @code python setup install --home=target_directory 
        @endcode to install the module to target_directory. This directory must be
        included in the @c PYTHONPATH environment variable 
 
Release notes
-------------

v1.0.0: initial stable version

