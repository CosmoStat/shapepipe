#!/usr/bin/env python

from distutils.core import setup

exec(open('shapepipe_base/info.py').read())

setup(
      name=__whoami__,
      packages=['shapepipe_base'],
      version=__version__,
      license="GNU GPLv3",
      description='Base utilties for ShapePipe modules.',
      author='Samuel Farrens',
      author_email='samuel.farrens@cea.fr',
      platforms=['posix', 'mac os'],
      classifiers=["Programming Language :: Python",
                   "Development Status :: 4 - Beta",
                   "License :: OSI Approved :: GNU General Public License "
                   "v3 (GPLv3)", "Operating System :: POSIX",
                   "Operating System :: MacOS", ],
     )
