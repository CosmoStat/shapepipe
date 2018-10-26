#! /usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

package_info = {}
infopath = os.path.abspath(os.path.join(os.path.dirname(__file__),
                           'shapepipe', 'info.py'))
with open(infopath) as open_file:
    exec(open_file.read(), package_info)

setup(
    name=package_info['__name__'],
    author=package_info['__author__'],
    author_email=package_info['__email__'],
    version=package_info['__version__'],
    packages=find_packages(),
    install_requires=['joblib',
                      'modopt',
                      'numpy'],
    description='Galaxy shape measurement pipeline.',
    long_description=package_info['__about__'],
    setup_requires=['pytest-runner', ],
    tests_require=['pytest', 'pytest-cov', ],
)
