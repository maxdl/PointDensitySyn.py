#!/usr/bin/env python3

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from os.path import join, dirname
from pointdensitysyn.version import version as __version__


setup(
    name="PointDensitySyn.py",
    version=__version__,
    description="Tool for analysis of immunogold labelling",
    long_description=open(join(dirname(__file__), "README.rst")).read(),
    author="Max Larsson",
    author_email="max.larsson@liu.se",
    license="MIT",
    url="http://www.liu.se/medfak/forskning/larsson-max/software",
    packages=find_packages(),
    entry_points={
    'console_scripts':
        ['PointDensitySyn = PointDensitySyn:main'],
    'gui_scripts':
        ['PointDensitySyn = PointDensitySyn:main']        
    },
    data_files=[('pointdensitysyn', ['pointdensitysyn/pds.ico'])],
    install_requires=['wxpython', 'openpyxl']
)