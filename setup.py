# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import sys

try:
    from setuptools import setup
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.", file=sys.stderr)


# load version info
exec(open("mageck_vispr/version.py").read())


setup(
    name="mageck-vispr",
    version=__version__,
    author="Johannes Köster",
    author_email="koester@jimmy.harvard.edu",
    description="MAGeCK-VISPR is a comprehensive quality control, analysis and "
    "visualization pipeline for CRISPR/Cas9 screens.",
    license="MIT",
    url="https://bitbucket.org/johanneskoester/mageck-vispr",
    packages=["mageck_vispr"],
    include_package_data=True,
    zip_safe=False,
    install_requires=["jinja2"],
    entry_points={"console_scripts": ["mageck-vispr = mageck_vispr.cli:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        #"Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
