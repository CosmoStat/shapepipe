#! /usr/bin/env python

from setuptools import setup, find_packages
import os

package_info = {}
infopath = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "shapepipe", "info.py")
)
with open(infopath) as open_file:
    exec(open_file.read(), package_info)


# Find scripts
def find_scripts():

    scripts = []
    scripts_dir = package_info["__scripts_dir__"]
    valid_extensions = package_info["__scripts_ext__"]

    sub_dirs = (
        os.path.join(scripts_dir, sub_dir)
        for sub_dir in os.listdir(scripts_dir)
        if os.path.isdir(os.path.join(scripts_dir, sub_dir))
    )

    for sub_dir in sub_dirs:
        scripts.extend(
            [
                os.path.join(sub_dir, val)
                for val in os.listdir(sub_dir)
                if os.path.isfile(os.path.join(sub_dir, val))
                and "__init__" not in val
                and val.endswith(valid_extensions)
            ]
        )

    return scripts


setup(
    name=package_info["__name__"],
    author=package_info["__author__"],
    author_email=package_info["__email__"],
    version=package_info["__version__"],
    packages=find_packages(),
    scripts=find_scripts(),
    setup_requires=package_info["__setups__"],
    install_requires=package_info["__installs__"],
    description="Galaxy shape measurement pipeline.",
    long_description=package_info["__about__"],
    tests_require=package_info["__tests__"],
)
