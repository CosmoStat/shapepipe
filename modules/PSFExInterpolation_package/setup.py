#!/usr/bin/env python

from distutils.core import setup

exec(open('PSFExInterpolation/info.py').read())

setup(
      name='PSFExInterpolation',
      packages=['PSFExInterpolation'],
      version=__version__,
      license="GNU GPLv3",
      description='Computation of PSFEx model at desired postion',
      author='Marc Gentile + Samuel Farrens',
      author_email='marc.gentile@epfl.ch',
      long_description=open('README.txt').read(),
      platforms=['posix', 'mac os'],
      classifiers=["Programming Language :: Python",
                   "Development Status :: 4 - Beta",
                   "License :: OSI Approved :: GNU General Public License "
                   "v3 (GPLv3)", "Operating System :: POSIX",
                   "Operating System :: MacOS", ],
     )
