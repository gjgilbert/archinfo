from __future__ import division, print_function

import sys

import setuptools
from setuptools import setup

sys.path.insert(0, "archinfo")
from version import __version__

long_description = \
"""
Archinfo is a python package to calculate measures describing the architectures of exoplanetary systems. The measures are described in [paper arXiv and ADS links]. Read the documentation at [website link].
"""

setup(
    name='archinfo'
    version=__version__,
    liscense='MIT'
    author='Gregory J. Gilbert'
    author_email='gjgilbert@uchicago.edu'
    packages=['archinfo']
    include_package_data=True
    url='http://github.com/gjgilbert/archinfo'
    description='Describing exoplanetary system architectures using information theory'
    long_description=long_description
    install_requires=['scipy', 'warnings']
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.0',]
)