#!/usr/bin/env python

from distutils.core import setup

exec(open('scatalog/scatalog_version.py').read())
setup(
      name='scatalog',
      version=__version__,
      license = "GNU GPLv3",
      description='Catalog Management',
      author='Marc Gentile',
      author_email='marc.gentile@epfl.ch',
      #url='http://www.python.org/sigs/sconfig-sig/',
      long_description=open('README.txt').read(),
      platforms=['posix','mac os'],
      packages=['scatalog']
     )


