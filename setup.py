# !/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""
setup.py for skysearcher.
"""

from numpy import get_include
from setuptools import setup

setup(
    name="skysearcher",
    version="0.0.1",
    description=(
        "Parallel group finding aslgorithm for finding stellar streams"
        "in the barionic fraction of simulated galactic dark matter halos"
    ),
    keywords='parallel stellar groupfinder',
    classifiers=[

        "Development Status :: 1 - Planning",
        # "Development Status :: 2 - Pre-Alpha",
        # "Development Status :: 3 - Alpha",
        # "Development Status :: 4 - Beta",
        # "Development Status :: 5 - Production/Stable",
        # "Development Status :: 6 - Mature",
        # "Development Status :: 7 - Inactive",

        # "Environment :: MacOS X",
        # "Environment :: Win32 (MS Windows)",
        "Environment :: X11 Applications",

        "Intended Audience :: Science/Research",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Other Audience",

        "Natural Language :: English",

        # "Operating System :: MacOS :: MacOS X",
        # "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",

        # "Programming Language :: Cython",
        # "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",

        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",

    ],
    author="Sol W Courtney",
    author_email="swc2124@columbia.edu",
    maintainer="Sol W Courtney",
    maintainer_email="swc2124@columbia.edu",
    url="https://github.com/swc2124/skysearcher",
    download_url="https://github.com/swc2124/skysearcher.git",
    license="MIT",
    packages=["skysearcher"],
    python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4",
    install_requires=[
        "numpy>=1.13.3",
        "mpi4py>=2.0.0",
        "numba>=0.36.2",
        "astropy>=2.0.3",
        "matplotlib>=2.1.1",
        "h5py>=1.10.1"
    ],
    include_package_data=True,
    include_dirs=[get_include()],
    package_data={
        "skysearcher": [
            "data/*.npy",
            "data/grids/*.npy",
            "data/plots/",
            "data/tables/*.hdf5",
            "data/tables/groupfinder/*.hdf5",
            "data/tables/groupfinder/mpi/*.hdf5"
        ]
    },
    entry_points={
        "console_scripts": [
            "skysearcher-newcfg = skysearcher.cmdlinetool:newcfg",
            "skysearcher-run = skysearcher.cmdlinetool:run",
        ],
    },

)
