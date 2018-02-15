#!/usr/bin/env python

from distutils.core import setup

setup(
      name='mkpsf',
      packages=['mkp', 'mkp.methods'],
      version='0.5.0',
      license = "GNU GPLv3",
      description='PSF reconstruction tool',
      author='Marc Gentile',
      author_email='marc.gentile@epfl.ch',
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


