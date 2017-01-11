#!/usr/bin/env python

from distutils.core import setup

exec(open('multifit/multifit_version.py').read())
setup(
      name='multifit',
      version=__version__,
      license = "GNU GPLv3",
      description='Multipurpose Fitting Frawework (MultiFit)',
      author='Marc Gentile',
      author_email='marc.gentile@epfl.ch',
      #url='http://www.python.org/sigs/sconfig-sig/',
      long_description=open('README.txt').read(),
      platforms=['posix','mac os'],
      packages=['multifit', 'multifit.modules', 'multifit.modules.models', 'multifit.modules.methods'],
     )


