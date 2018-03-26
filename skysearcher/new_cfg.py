# !/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


"""
===========================
New onfiguration file maker
===========================

This script contains a single function ``new_rc()``.  It is to be called
in the event that no configuration file exists in the main skysearcher
directory ``skysearcher/skysearcher``.  After running this function, there
will be a new configuration file titled ``"rc.cfg"``.  Skysearcher will read
this file at the start of each run.

note: The file can be manually edited as text prior to running
skysearcher.
"""

from warnings import warn
import configparser
import os


__all__ = ["new_rc"]


def new_rc(rc_fh=None, output=True, tbl_ext="hdf5"):
    """Create a new rc configuration file to the provided file handle.

    Parameters
    ----------
    rc_fh : None, optional
        File name of the output config file.
    tbl_ext : str, optional
        Description
    output : bool, optional
        If True then print newly generated configfile contents to sdout.

    Keyword Arguments
    -----------------
    rc_fh {str} -- full path and name of rc file (default: {None})
    tbl_ext {str} -- extention for tables. (default: {"hdf5"})

    Returns
    -------
    rc file handle str.

    Example
    -------
        >>> import ConfigParser
        >>> config = ConfigParser.RawConfigParser()
        >>> config.read("rc.cfg")

    .. code-block:: cfg

        user@machine$ python new_cfg.py
        [ PATH ]
        data_dir = ~/$USER/skysearcher/data
        plot_dir = ~/$USER/skysearcher/data/plots
        table_dir = ~/$USER/skysearcher/data/tables
        mpi_table_dir = ~/$USER/skysearcher/data/tables/groupfinder/mpi
        grid_dir = ~/$USER/skysearcher/data/grids
        grid_file_designator = grid
        grid_ext = npy
        table_file_designator = table
        table_format = hdf5
        table_hdf5_path = data
        table_ext = .hdf5

        [ Search Extent ]
        r_start = 5
        r_stop = 285
        r_step = 1
        annulus_scale = 0.05
        annulus_phi_step = 720

        [ Accept Reject ]
        xbox_cut = 0.1
        min_log_nstars = 2.00
        min_n_segments = 4
        n_skips = 2

        [ Run Time ]
        save_interval = 2

        [ Data ]
        d_mpc = 4.0


    .. seealso::
    
        ConfigParser: http://docs.python.org/2/library/configparser.html
    """
    
    # Make a new configuration file parser object (config).
    config = configparser.ConfigParser()
    
    # Get a string value for:
    # 1.) path to parent directory (pardir_path).
    # 2.) name of current directory (curdir_name).    
    if not rc_fh:

        # Make a new file name (file handle or fh).
        rc_fh = "rc"

    rc_fh += os.path.extsep + "cfg"

    try:
        import skysearcher as ss
    except ImportError as err:
        warn("\n--> Skysearcher is not installed!\n")
        # Get a string value for the path to the
        # current directory (curdir_path).
        data_path = os.path.abspath(os.path.curdir)
        rc_fh_path = os.path.join(data_path, rc_fh)
    else:
        data_path = os.path.dirname(ss.__file__)
        rc_fh_path = os.path.join(data_path, rc_fh)


    msg = "saving [ " + rc_fh + " ] to : " + rc_fh_path
    print("\n[ new_rc ] [ Configuration file location ]")
    print("-"*len(msg))
    print(msg)
    print("-"*len(msg), "\n")

    # Name the directory PATH for all data folders (data_dir).
    data_dir = os.path.join(data_path, "data")
    table_dir = os.path.join(data_dir, "tables")
    table_format = tbl_ext

    config["PATH"] = {
        "data_dir": data_dir,
        "plot_dir": os.path.join(data_dir, "plots"),
        "table_dir": table_dir,
        "mpi_table_dir": os.path.join(table_dir, "groupfinder", "mpi"),
        "grid_dir": os.path.join(data_dir, "grids"),
        "grid_file_designator": "grid",
        "grid_ext": "npy",
        "table_file_designator": "table",
        "table_format": table_format,
        "table_hdf5_path": "data",
        "table_ext": table_format}

    config["Search Extent"] = {
        "r_start": "1",
        "r_stop": "285",
        "r_step": "1",
        "annulus_scale": "0.05",
        "annulus_phi_step": "720"}

    config["Accept Reject"] = {
        "xbox_cut": "0.1",
        "min_log_nstars": "2.25",
        "min_n_segments": "4",
        "n_skips": "2"}

    config["Run Time"] = {
        "save_interval": "2"}

    config["Data"] = {
        "d_mpc": "4.0"}
     
    # Open file object (configfile) to write to.
    with open(rc_fh_path, "w") as configfile:
        config.write(configfile)

    if output:
        msg = str("[ "+ rc_fh_path + " ]")
        print("-"*len(msg))
        print(msg)
        print("-"*len(msg))
        for sec in config.sections():
            print("[", sec, "]")
            for opt in config.options(sec):
                print(opt, "=", config.get(sec, opt))
            print(" ")
    return rc_fh


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        prog="new_conf.py",
        description="MPI executable function")
    parser.add_argument("-f", "--filename",
        type=str,
        default="rc",
        choices=["rc", "conf", "configfile"],
        help="name of configuration file without extention")
    parser.add_argument("-o", "--output",
        type=bool,
        default=True,
        choices=[True, False],
        help="increase output verbosity")
    parser.add_argument("-tbl-ext", "--table-extention",
        type=str,
        default="hdf5",
        choices=["hdf5"],
        help="extention for tables")
    args = parser.parse_args(os.sys.argv[1:])

    new_rc(rc_fh=args.filename,
        tbl_ext=args.table_extention,
        output=args.output)
