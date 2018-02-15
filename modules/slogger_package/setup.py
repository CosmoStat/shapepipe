#!/usr/bin/env python

from distutils.core import setup

exec(open('slogger/version.py').read())
setup(
      name='slogger',
      version=__version__,
      license = "GNU GPLv3",
      description='Lightweight File Logging facility',
      author='Marc Gentile',
      author_email='marc.gentile@epfl.ch',
      #url='http://www.python.org/sigs/sconfig-sig/',
      long_description=open('README.txt').read(),
      platforms=['posix','mac os'],
      packages=['slogger']
     )


