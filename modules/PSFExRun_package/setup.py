#!/usr/bin/env python

# !!! Don't forget to update the package name!

from distutils.core import setup

exec(open('PSFExRun/info.py').read())

setup(
      name='PSFExRun',
      packages=['PSFExRun'],
      version=__version__,
      license="GNU GPLv3",
      description='Parallel execution of multiple PSFEx processes',
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
