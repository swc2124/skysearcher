# !/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

"""
Summary
"""

import os
import argparse

import skysearcher as ss


def newcfg():
    """
    Summary
    """
    from .new_cfg import new_rc

    parser = argparse.ArgumentParser(prog="new_conf.py",
         description="MPI executable function")
    parser.add_argument("-f", "--filename",
                        type=str,
                        default="rc",
                        help="name of configuration file without extention")
    parser.add_argument("-o", "--output",
                        type=bool,
                        default=True,
                        help="increase output verbosity")
    args = parser.parse_args(os.sys.argv[1:])

    new_rc(rc_fh=args.filename,
           output=args.output)


def run():
    """
    Summary
    """
    from subprocess import call

    print(ss.__doc__)

    parser = argparse.ArgumentParser(
        prog="mpi_search.py",
        description="""MPI executable function. \
        mpiexec -n python mpi_search.py""")
    parser.add_argument("-n", "--nproc",
                        type=int,
                        default=1,
                        help="number of parallel processes to run")
    _args = parser.parse_args(os.sys.argv[1:])
    _path = os.path.join(os.path.dirname(ss.__file__), parser.prog)
    try:
        call(args=["mpiexec", "-n", str(_args.nproc), "python", _path])

    except KeyboardInterrupt as err:
        print("\n -- \n")
        print(err)
        os.sys.exit(0)

