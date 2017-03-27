#!/usr/bin/env python

from distutils.core import setup

exec(open('sconfig/sconfig_version.py').read())
setup(
      name='sconfig',
      version=__version__,
      license = "GNU GPLv3",
      description='Configuration file parser',
      author='Marc Gentile',
      author_email='marc.gentile@epfl.ch',
      #url='http://www.python.org/sigs/sconfig-sig/',
      long_description=open('README.txt').read(),
      platforms=['posix','mac os'],
      classifiers = [
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX",
        "Operating System :: Mac"
      ],
      packages=['sconfig']
     )


