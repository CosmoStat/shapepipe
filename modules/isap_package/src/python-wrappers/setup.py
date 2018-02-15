#!/usr/bin/env python

from distutils.core import setup

setup(
      name='IsapPyWrapper',
      packages=['IsapPyWrapper','IsapPyWrapper.sparse2d'],
      version='0.1',
      license = "GNU GPLv3",
      description='Wrappers for ISAP GREAT3',
      author='Florent Sureau',
      author_email='florent.sureau@cea.fr',
      #url='http://www.python.org/sigs/quadg3-sig/',
      long_description=open('README.txt').read(),
      platforms=['posix','mac os'],
      classifiers = [
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX",
        "Operating System :: MacOS", 
        ],
     )


