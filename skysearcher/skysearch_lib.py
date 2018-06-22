# !/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""
Attributes
----------
PLOT_DIR : str
    Path to plot output directory.

TABLE_DIR : str
    Path to hdf5 table input and output directory.

MPI_TABLE_DIR : str
    Path to mpi hdf5 table directory within TABLE_DIR.

GRID_DIR : str
    Path to numpy array input files directory.

TABLE_FORMAT : str
    Table format for astropy.tables.Table (default: {"hdf5"})

TABLE_HDF5_PATH : str
    Table path for hdf5 table (default: {"hdf5"})

TABLE_EXT : str
    Table file-type extension (default: {"hdf5"})

GRID_EXT : str
    Extension for array files. (default: {"npy"})

R_START : int
    Starting radius in units of Kpc to start search from.

R_STOP : int
    Ending search radius in units of Kpc.

R_STEP : int
    Step interval in units of Kpc

R_SCALE : float
    Percent to multiply radius by to get + and - Kpc for annulus.

ANNULUS_PHI_STEP : float
    How many sections to divide each annulus by.

XBOX_CUT : float
    Minimum contrast density a feature must have.

MIN_LOG_NSTARS : float
    Minimum log10 value for number of stars a feature must have.

MIN_N_SEGMENTS : int
    Minimum number of segments a feature must be

N_SKIPS : int
    Allowed number of skipped segments before ending feature.

SAVE_INTERVAL : int
    Number of features to hold between saves.

D_MPC : float
    Distance to target halo in units of Mpc.
"""

import configparser
import os
import copy

from time import sleep
from time import time
import datetime

import numpy as np

from mpi4py import MPI

from numba import jit

from astropy.table import Table
# from astropy.table import vstack

__all__ = [
    "np", "os", "time", "clear_tables", "pause", "clear_tables",
    "grid_list", "kpc_to_arcmin", "record_table", "radii", "TABLE_FORMAT",
    "save_record_table", "fix_rslice", "load_grid", "satid_setup",
    "mu_idx", "mu_idx2", "get_annuli", "get_idx", "get_xbox", "COMM",
    "TABLE_HDF5_PATH", "new_sat_stars", "count_strs", "dom_satid",
    "MPI_PROC_NAME", "MPI_SIZE", "MPI_RANK", "STDOUT", "D_MPC", "TABLE_EXT",
    "TABLE_COLUMNS", "DATA_DIR", "TABLE_DIR", "PLOT_DIR", "MIN_LOG_NSTARS",
    "START_TIME", "GRID_DIR", "R_START", "R_STOP", "R_STEP", "N_SKIPS",
    "MIN_N_SEGMENTS", "SAVE_INTERVAL", "XBOX_CUT", "MPI_TABLE_DIR",
    "R_SCALE", "ANNULUS_PHI_STEP", "report", "datetime"
]

COMM = MPI.COMM_WORLD
MPI_RANK = COMM.Get_rank()
MPI_SIZE = COMM.Get_size()
MPI_PROC_NAME = MPI.Get_processor_name()
START_TIME = time()
STDOUT = os.sys.stdout

def sortout_rcfile(name="rc", ext="cfg"):
    """
    Return the correct file name for the skysearcher rc.cfg file.

    Make sure we"re in skysearcher/skysearcher & save path variables.
    Determine if there exists a configuration file by making a list
    of filenames for all the files in the current directory. Isolate all
    files with .cfg extension within curdir_files.  Ask if there is at least
    one configuration file present?
    - If there is:
    Ask if there is more than one configuration file.  If there is then
    print out available options and ask user for selection.
    - If there is not:
    Use the configuration file generator function new_cfg.new_rc().

    Parameters
    ----------
    name : str, optional
        config file name
    ext : str, optional
        config file extension to be used

    Returns
    -------
    str
        the proper file name for the skysearcher config file

    """

    # configuration file settings
    rcfile_ext = os.path.extsep + ext
    new_rc_fh = name + rcfile_ext

    try:
        import skysearcher as ss
    except ImportError as err:
        from warnings import warn
        warn("\n --> Skysearcher is not installed!\n")
        rc_fh_path = os.path.abspath(os.path.curdir)
    else:
        rc_fh_path = os.path.dirname(ss.__file__)

    msg = "Using " + rc_fh_path + " as top level dir"
    print("\n[ sortout_rcfile ] [ Configuration file location ]")
    print("-" * len(msg))
    print(msg)
    print("-" * len(msg), "\n")

    # Determine if there exists a configuration file by making a list
    # of filenames for all the files in the current directory
    # (curdir_files).
    curdir_files = os.listdir(rc_fh_path)

    # Isolate all files with .cfg extension within curdir_files
    # (rc_files).
    rc_files = [f for f in curdir_files
                if os.path.splitext(f)[1] == rcfile_ext]

    # Ask if there is at least one configuration file present?
    # If there is:
    if len(rc_files):

        # Ask if there is more than one configuration file.
        # If there is:
        if len(rc_files) > 1:

            # Print out available options.
            for i, rc_filename in enumerate(rc_files):
                print("     [", i, "]   ", rc_filename)

            # Ask user for selection.
            selection = int(input("\n\n       Enter selection: "))
            rc_fh = rc_files.pop(selection)

        # If there is only one:
        else:
            # Take the only name in the list.
            rc_fh = rc_files.pop()

        rc_fh = os.path.join(rc_fh_path, rc_fh)

    # If there is not:
    else:
        # Use the configuration file generator function new_cfg.new_rc()
        # and clean up.
        try:
            import new_cfg
        except ImportError:
            from . import new_cfg
        print("[ sortout_rcfile ] [ Making new rc.cfg ]")
        rc_fh = new_cfg.new_rc(rc_fh=new_rc_fh)

    print(os.path.abspath(rc_fh))

    return rc_fh
# <configuration file values>
# Make a new configuration file parser object (config).
_CONFIG = configparser.ConfigParser()

# Read in rc.file and unpack all values.
_CONFIG.read(sortout_rcfile())
DATA_DIR = _CONFIG.get("PATH", "data_dir")
TABLE_DIR = _CONFIG.get("PATH", "table_dir")
PLOT_DIR = _CONFIG.get("PATH", "plot_dir")
GRID_DIR = _CONFIG.get("PATH", "grid_dir")
MPI_TABLE_DIR = _CONFIG.get("PATH", "mpi_table_dir")
GRID_EXT = _CONFIG.get("PATH", "grid_ext")
TABLE_EXT = _CONFIG.get("PATH", "table_ext")
TABLE_FORMAT = _CONFIG.get("PATH", "table_format")
TABLE_HDF5_PATH = _CONFIG.get("PATH", "table_hdf5_path")
R_START = _CONFIG.getint("Search Extent", "r_start")
R_STOP = _CONFIG.getint("Search Extent", "r_stop")
R_STEP = _CONFIG.getint("Search Extent", "r_step")
R_SCALE = _CONFIG.getfloat("Search Extent", "annulus_scale")
ANNULUS_PHI_STEP = _CONFIG.getint("Search Extent", "annulus_phi_step")
_ANNULUS_PHI_STEP = None
XBOX_CUT = _CONFIG.getfloat("Accept Reject", "xbox_cut")
MIN_LOG_NSTARS = _CONFIG.getfloat("Accept Reject", "min_log_nstars")
MIN_N_SEGMENTS = _CONFIG.getint("Accept Reject", "min_n_segments")
N_SKIPS = _CONFIG.getint("Accept Reject", "n_skips")
SAVE_INTERVAL = _CONFIG.getint("Run Time", "save_interval")
D_MPC = _CONFIG.getfloat("Data", "d_mpc")
TABLE_COLUMNS = [

    # Halo.
    ("halo", "i"),

    # Annulus location values.
    ("radius", "i"),
    ("r0", "f"),
    ("r1", "f"),
    ("annuli_step", "f"),
    # ("n_empty_segments", "i"),

    # Annulus content values.
    # ("log_n_boxes_in_ann", "i"),
    # ("log_n_stars_in_ann", "i"),
    ("Log10(mu)", "f"),

    # Feature content values.
    ("xbox_max", "f"),
    ("Log10(n_stars_max)", "f"),
    ("domsat_purity", "f"),
    ("domsat_id", "i"),
    ("domsat_sig", "f"),
    ("nsats", "f"),
    ("domsat_mass", "f"),
    ("domsat_atime", "f"),
    ("domsat_j", "f"),

    # Feature location values.
    ("deg0", "f"),
    # ("deg1", "f"),
    ("extent", "f"),
    # ("n_segments", "i"),
    # ("n_boxes", "i"),

    # ("MPI_RANK", "i")
]


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)


def ct(text, colour=WHITE):
    return "\x1b[1;%dm" % (30 + colour) + text + "\x1b[0m"


def render_and_pad(reqd_width, components, sep=""):
    temp = []
    actual_width = 0
    for fmt_code, text in components:
        actual_width += len(text)
        strg = "\x1b[1;%dm" % (30 + fmt_code) + text + "\x1b[m"
        temp.append(strg)
    if temp:
        actual_width += len(temp) - 1
    npad = reqd_width - actual_width
    if npad < 0:
        npad = 0
    return sep.join(temp) + " " * npad


def clear_tables(_target_dir):
    """
    Clear all MPI temp tables from :term:`MPI_TABLE_DIR`
    before creating new ones.

    Parameters
    ----------
    _target_dir : str, optional
        PATH to target directory.

    Raises
    ------
    IOError
        In case they cant be deleted.
    TypeError
        :term:`MPI_TABLE_DIR` may not be defined.

    Example
    -------
        >>> from skysearcher.skysearch_lib import *
        >>> import os
        >>> os.listdir(MPI_TABLE_DIR)
        ['17.hdf5',
         '2.hdf5',
         '11.hdf5',
         ...,
         '0.hdf5',
         '15.hdf5',
         '12.hdf5']
        >>> clear_tables()
        >>> os.listdir(MPI_TABLE_DIR)
        []
    """
    try:
        for fh in os.listdir(_target_dir):
            os.remove(os.path.join(_target_dir, fh))
    except FileNotFoundError as err:
        _CONFIG = configparser.ConfigParser()
        _CONFIG.read(sortout_rcfile())
        for itm in _CONFIG.items("PATH"):
            if itm[0].endswith("_dir"):
                if not os.path.isdir(itm[1]):
                    os.mkdir(itm[1])
    except IOError as err:
        print(err)
    except TypeError as err:
        print(err)


def pause(rnk=MPI_RANK):
    """
    Slow down MPI process according to their :term:`rank`.

    Parameters
    ----------
    rnk : int, optional
        Integer representing MPI process rank

    Example
    -------
        >>> skysearcher.skysearch_lib import pause
        >>> pause(rnk=1)
    """
    sleep(MPI_RANK * 0.025)


def grid_list():
    """
    Load file handles for numpy data arrays.

    Provide a list of all available numpy data arrays (:term:`grids`).

    Returns
    -------
    list
        List of grids to be used in search.
        (halo,  PATH)

    Example
    -------
        >>> from skysearcher.skysearch_lib import *
        >>> grids = grid_list()
        >>> grids
        [('halo07', '$PATH/skysearcher/data/grids/halo07_4.0Mpc_h158_grid.npy'),
         ('halo09', '$PATH/skysearcher/data/grids/halo09_4.0Mpc_h158_grid.npy'),
         ('halo20', '$PATH/skysearcher/data/grids/halo20_4.0Mpc_h158_grid.npy'),
         ...
         ('halo08', '$PATH/skysearcher/data/grids/halo08_4.0Mpc_h158_grid.npy'),
         ('halo15', '$PATH/skysearcher/data/grids/halo15_4.0Mpc_h158_grid.npy'),
         ('halo02', '$PATH/skysearcher/data/grids/halo02_4.0Mpc_h158_grid.npy')]
    """
    return [
        (fh.split("_")[0], os.path.join(GRID_DIR, fh))
        for fh in os.listdir(GRID_DIR)
        if fh.endswith(GRID_EXT)]


@jit(cache=True)
def kpc_to_arcmin(d_mpc=D_MPC):
    """
    Make the :term:`mod` for the distance.

    Convert :term:`Kpc` to :term:`Arc-min`.

    Parameters
    ----------
    d_mpc : float, optional
        Distance in :term:`Mpc`.

    Returns
    -------
    float
        Ratio of :term:`Kpc` to :term:`Arc-min`.

    Example
    -------
        >>> from skysearcher.skysearch_lib import kpc_to_arcmin
        >>> mod = kpc_to_arcmin(d_mpc=4.0)
        >>> mod
        0.7386310975322501
    """
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    arcmin_mod = 3437.746  # (60.0 * 180.0) / np.pi
    return np.square(arcmin_mod * (d / D))


def record_table(_names=TABLE_COLUMNS, _meta=True):
    """
    Load the record table (:term:`r_table`).

    Make a new astropy table to use for MPI output.

    Parameters
    ----------
    _names : list, optional
        List of column names and data types to be used.
        (column name, dtype)
    _meta : bool, optional
        Attach meta data to table or not.

    Returns
    -------
    astropy.tables.Table
        Record table to use for storing MPI data.

    Examples
    --------
    >>> from skysearcher.skysearch_lib import *
    >>> record_table(_names=TABLE_COLUMNS)
    <Table length=0>
     halo radius    r0      r1   annuli_step ... domsat_atime domsat_j   deg0   extent
    int32 int32  float32 float32   float32   ...   float32    float32  float32 float32
    ----- ------ ------- ------- ----------- ... ------------ -------- ------- -------

    """
    output_colums_ = []
    output_dtyps_ = []
    for name_, dtyp_ in _names:
        output_colums_.append(name_)
        output_dtyps_.append(dtyp_)
    r_tbl_ = Table(names=output_colums_, dtype=output_dtyps_)
    if _meta:
        r_tbl_.meta["r_start"] = R_START
        r_tbl_.meta["r_stop"] = R_STOP
        r_tbl_.meta["r_step"] = R_STEP
        r_tbl_.meta["r_scale"] = R_SCALE
        r_tbl_.meta["fh"] = os.path.join(
            MPI_TABLE_DIR,
            str(MPI_RANK) + os.path.extsep + TABLE_EXT)
    return r_tbl_


def radii():
    """
    Load list of radii (:term:`_radii`).
    i.e radii[x] = (:term:`r`, :term:`r_start`, :term:`r_stop`)

    The list of radii to search.  Set from rc.cfg file.

    Returns
    -------
    list
        List of radii in the form (center, inner, outer)

    Examples
    --------
    >>> from skysearcher.skysearch_lib import radii
    >>> _radii = radii()
    >>> _radii
    [(5, 4.75, 5.25),
     (6, 5.7, 6.3),
     (7, 6.65, 7.35),
     ...
     (11, 10.45, 11.55),
     (12, 11.4, 12.6),
     (13, 12.35, 13.65)]

    """
    radii_list_ = []
    for r in list(range(R_START, R_STOP, R_STEP)):
        annulus_scale_ = R_SCALE * r
        annulus_scale_ *= (1.25 - (r / R_STOP)) # 0.6
        radii_list_.append((r, r - annulus_scale_, r + annulus_scale_))
        """
        start = r * np.tan(45)
        stop = (r + R_STEP) * np.tan(45)
        radii_list_.append((start, start, stop))
        """
    return radii_list_


def save_record_table(_table, _rnk=MPI_RANK):
    """
    Save record hdf5 table.

    Parameters
    ----------
    _table : astropy.table.Table
        Data table for storing search data
    _rnk : int, optional
        :term:`MPI_RANK`

    Raises
    ------
    IOError
        Description
    TypeError
        Missing 1 required positional argument: '_table'

    Examples
    --------
    >>> from skysearcher.skysearch_lib import *
    >>> r_table = record_table(_names=TABLE_COLUMNS)
    >>> save_record_table(_table=r_table)

    """
    try:
        _table.write(_table.meta["fh"],
                     format=TABLE_FORMAT,
                     path=TABLE_HDF5_PATH,
                     overwrite=True,
                     append=True)
    except IOError:
        pass
        # raise IOError()
    except TypeError:
        pass
        # raise TypeError()
    except OSError:
        print(ct("[ MPI_RANK : " + str(MPI_RANK) + " FAILED TO SAVE ]", RED))


def fix_rslice(_grid, _dtype=np.float32):
    """
    Summary

    Parameters
    ----------
    _grid : TYPE
        Description
    _dtype : TYPE, optional
        Description

    Returns
    -------
    TYPE
        Description

    Raises
    ------
    IndexError
        Description

    """
    grid_ = np.concatenate(
        (_grid.astype(_dtype), np.zeros(
            (600, 600, 1), dtype=_dtype)), axis=2)
    center = 300.0 - 0.5
    try:
        for i in range(grid_.shape[0]):
            for q in range(grid_.shape[1]):
                value = np.sqrt(
                    np.square(i - center, dtype=_dtype) +
                    np.square(q - center, dtype=_dtype),
                    dtype=_dtype)
                grid_[i, q, 5] = value
                grid_[i, q, 6] = np.arctan2(
                    q - center, i - center, dtype=_dtype)
    except IndexError as err:
        raise IndexError(err)
    else:
        return grid_


def load_grid(_grid_fh):
    """
    Get the data :term:`grid` for this halo (:term:`grid`).

    Parameters
    ----------
    _grid_fh : str
        File handle (absolute PATH) for a grid file.

    Returns
    -------
    list
        List of Tuples [(halo_name, halo_data_path), ...]

    Examples
    --------
    >>> from skysearcher.skysearch_lib import *
    >>> grids = grid_list()
    >>> grid_info = grids[0]
    >>> grid_info
    ('halo07',
     '$PATH/skysearcher/data/grids/halo07_4.0Mpc_h158_grid.npy')
    >>> halo, grid_fh = grid_info
    >>> halo
    'halo07'
    >>> grid_fh
    '$PATH/skysearcher/data/grids/halo07_4.0Mpc_h158_grid.npy'
    >>> grid = load_grid(grid_fh)
    >>> grid
    array([[[   0., 0., 0., ..., 0., 423.55697632, -2.3561945 ],
            [   0., 0., 0., ..., 0., 422.85043335, -2.35786676],
            [   0., 0., 0., ..., 0., 422.14511108, -2.35954452],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 2.35954452],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35786676],
            [   0., 0., 0., ..., 0., 423.55697632, 2.3561945 ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -2.35452223],
            [   0., 0., 0., ..., 0., 422.14276123, -2.3561945 ],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35787249],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 2.35787249],
            [   0., 0., 0., ..., 0., 422.14276123, 2.3561945 ],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35452223]
            ],
           [[   0., 0., 0., ..., 0., 422.14511108, -2.35284448],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35451674],
            [   0., 0., 0., ..., 0., 420.72854614, -2.3561945 ],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 2.3561945 ],
            [   0., 0., 0., ..., 0., 421.43624878, 2.35451674],
            [   0., 0., 0., ..., 0., 422.14511108, 2.35284448]],
            ...,
           [[   0., 0., 0., ..., 0., 422.14511108, -0.7887482 ],
            [   0., 0., 0., ..., 0., 421.43624878, -0.787076  ],
            [   0., 0., 0., ..., 0., 420.72854614, -0.78539819],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 0.78539819],
            [   0., 0., 0., ..., 0., 421.43624878, 0.787076  ],
            [   0., 0., 0., ..., 0., 422.14511108, 0.7887482 ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -0.78707045],
            [   0., 0., 0., ..., 0., 422.14276123, -0.78539819],
            [   0., 0., 0., ..., 0., 421.43624878, -0.78372031],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 0.78372031],
            [   0., 0., 0., ..., 0., 422.14276123, 0.78539819],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78707045]
            ],
           [[   0., 0., 0., ..., 0., 423.55697632, -0.78539819],
            [   0., 0., 0., ..., 0., 422.85043335, -0.78372592],
            [   0., 0., 0., ..., 0., 422.14511108, -0.78204811],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 0.78204811],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78372592],
            [   0., 0., 0., ..., 0., 423.55697632, 0.78539819]
           ]], dtype=float32)
    >>> grid.shape
    (600, 600, 7)

    """
    return fix_rslice(np.load(_grid_fh))


def satid_setup(_halo, attempts=0):
    """
    Get a list of :term:`satid`'s and a table for counting satids
    per region (:term:`satid_list`) (:term:`satid_table`).

    Parameters
    ----------
    _halo : str
        Halo name i.e. "halo02"
    attempts : int, optional
        No longer used.

    Returns
    -------
    tuple
        (list, astropy.table.Table)

    Raises
    ------
    OSError
        Unable to open file.  No such file or directory.
    TypeError
        Missing 1 required positional argument: '_halo'

    Examples
    --------

    >>> from skysearcher.skysearch_lib import *
    >>> halo = "halo07"
    >>> satid_list, satid_table = satid_setup(halo)
    >>> satid_list
    [215,
     217,
     218,
     ...,
     316,
     317,
     318]
    >>> satid_table
    <Table length=8243050>
    satids    Phi     Rads
    int16     rad     kpc
    uint16  float16  uint16
    ------ --------- ------
       225   -3.0566    229
       225   -3.0605    229
       254   -3.0664    229
       ...       ...    ...
       225  0.076111    229
       225  0.080444    229
       225  0.084778    229
    >>> satid_table.keys()
    ['satids', 'Phi', 'Rads']
    >>> satid_table["satids"].min() >= min(satid_list)
    True
    >>> satid_table["satids"].max() <= max(satid_list)
    True

    """
    table_fh = os.path.join(TABLE_DIR, _halo + os.path.extsep + TABLE_EXT)
    try:
        _tbl = Table.read(
            table_fh,
            format=TABLE_FORMAT,
            path=TABLE_HDF5_PATH)
    except TypeError as err:
        raise TypeError
    except OSError as err:
        raise OSError
    else:
        keeps = ["satids", "Rads", "Phi"]
        _tbl.keep_columns(keeps)
        return _tbl.meta["satids"], _tbl


def mu_idx(_grid, _r0, _r1):
    """
    Get the mean number of stars for this annulus (:term:`mu`) and an array of
    indices's representing the corresponding array segments (:term:`r_idx`).

    Parameters
    ----------
    _grid : numpy.ndarray
        The data grid for this halo (:term:`grid`)
    _r0 : float
        Starting radius (always < :term:`r`)
    _r1 : float
        Ending radius (always > :term:`r`)

    Returns
    -------
    tuple
        (float, numpy.ndarray)

    Examples
    --------
    >>> from skysearcher.skysearch_lib import *
    >>> grids = grid_list()
    >>> grid_info = grids[0]
    >>> halo, grid_fh = grid_info
    >>> grid = load_grid(grid_fh)
    >>> _radii = radii()
    >>> r, r_start, r_stop = _radii[0]
    >>> mu, r_idx = mu_idx(grid, r_start, r_stop)
    >>> mu
    6650.8335
    >>> r_idx
    array([[False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           ...,
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False]], dtype=bool)

    """
    grid_idx_ = np.logical_and(
        _grid[:, :, 5] >= _r0,
        _grid[:, :, 5] < _r1)
    _idx_ = np.nonzero(grid_idx_)
    if len(_idx_[0]):
        mean = _grid[:, :, 1][_idx_].mean()
    else:
        mean = 1.0
    return mean, grid_idx_


def mu_idx2(_grid, _r_indx, d0, d1):
    """
    Get the mean number of stars for this sub-annulus-section (:term:`mu`) and
    the array indices's for the corresponding array elements (:term:`r_idx2`).

    Parameters
    ----------
    _grid : numpy.ndarray
        The data grid for this halo (:term:`grid`)
    _r_indx : numpy.ndarray
        Indices's of array elements within the current :term:`annulus`.
    d0 : float
        Starting point (-pi, pi].
    d1 : float
        Ending point (-pi, pi].

    Returns
    -------
    tuple
        float, numpy.ndarray

    Example
    -------
    >>> from skysearcher.skysearch_lib import *
    >>> grids = grid_list()
    >>> grid_info = grids[0]
    >>> halo, grid_fh = grid_info
    >>> grid = load_grid(grid_fh)
    >>> grid
    array([[[   0., 0., 0., ..., 0., 423.55697632, -2.3561945 ],
            [   0., 0., 0., ..., 0., 422.85043335, -2.35786676],
            [   0., 0., 0., ..., 0., 422.14511108, -2.35954452],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 2.35954452 ],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35786676 ],
            [   0., 0., 0., ..., 0., 423.55697632, 2.3561945  ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -2.35452223],
            [   0., 0., 0., ..., 0., 422.14276123, -2.3561945 ],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35787249],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 2.35787249 ],
            [   0., 0., 0., ..., 0., 422.14276123, 2.3561945  ],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35452223 ]
            ],
           [[   0., 0., 0., ..., 0., 422.14511108, -2.35284448],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35451674],
            [   0., 0., 0., ..., 0., 420.72854614, -2.3561945 ],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 2.3561945  ],
            [   0., 0., 0., ..., 0., 421.43624878, 2.35451674 ],
            [   0., 0., 0., ..., 0., 422.14511108, 2.35284448 ]],
            ...,
           [[   0., 0., 0., ..., 0., 422.14511108, -0.7887482 ],
            [   0., 0., 0., ..., 0., 421.43624878, -0.787076  ],
            [   0., 0., 0., ..., 0., 420.72854614, -0.78539819],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 0.78539819 ],
            [   0., 0., 0., ..., 0., 421.43624878, 0.787076   ],
            [   0., 0., 0., ..., 0., 422.14511108, 0.7887482  ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -0.78707045],
            [   0., 0., 0., ..., 0., 422.14276123, -0.78539819],
            [   0., 0., 0., ..., 0., 421.43624878, -0.78372031],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 0.78372031 ],
            [   0., 0., 0., ..., 0., 422.14276123, 0.78539819 ],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78707045 ]
            ],
           [[   0., 0., 0., ..., 0., 423.55697632, -0.78539819],
            [   0., 0., 0., ..., 0., 422.85043335, -0.78372592],
            [   0., 0., 0., ..., 0., 422.14511108, -0.78204811],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 0.78204811 ],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78372592 ],
            [   0., 0., 0., ..., 0., 423.55697632, 0.78539819 ]
           ]], dtype=float32)
    >>> grid.shape
    (600, 600, 7)
    >>> _radii = radii()
    >>> r, r_start, r_stop = _radii[0]
    >>> r
    5
    >>> r_start
    4.75
    >>> r_stop
    5.25
    >>> mu, r_idx = mu_idx(grid, r_start, r_stop)
    >>> mu
    6650.8335
    >>> r_idx
    array([[False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           ...,
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False]], dtype=bool)
    >>> annuli, annuli_step = get_annuli(r)
    >>> annuli
    [(-3.1730085801256909, -3.1642795331377247),
     (-3.1642795331377247, -3.1555504861497585),
     (-3.1555504861497585, -3.1468214391617924),
     ...,
     (3.1468214391617924, 3.1555504861497585),
     (3.1555504861497585, 3.1642795331377247),
     (3.1642795331377247, 3.1730085801256909)]
    >>> annuli_step
    0.0087290469879661367
    >>> _deg0, _deg1 = annuli[0]
    >>> _deg0
    -3.1730085801256909
    >>> _deg1
    -3.1642795331377247
    >>> mu, r_idx2 = mu_idx2(grid, r_idx, _deg0, _deg1)
    >>> mu
    1.0
    >>> r_idx2
    array([[False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           ...,
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False]], dtype=bool)

    """
    # TODO add this to rc.cfg file.
    space = np.pi / 5.0
    _d0 = d0 - space
    _d1 = d1 + space
    if _d0 < -np.pi:
        _d0 += (2.0 * np.pi)
    if _d1 > np.pi:
        _d1 -= (2.0 * np.pi)
    grid_idx_ = np.logical_and(
        np.logical_and(
            _grid[:, :, 6] >= _d0,
            _grid[:, :, 6] < _d1),
        _r_indx)
    _idx_ = np.nonzero(grid_idx_)
    if len(_idx_[0]):
        mean = _grid[:, :, 1][_idx_].mean()
    else:
        mean = 1.0
    return mean, grid_idx_


def get_annuli(_r):
    """
    Load array of annuli segments and step value (:term:`annuli`) &
    (:term:`annuli_step`).

    Parameters
    ----------
    _r : int
        Radius Kpc.

    Returns
    -------
    tuple
        (list of tuples, float)
            list of tuples:
            annuli_step = float
            annuli[x] = (deg_0, deg_1)
            annuli, annuli_step

    Example
    -------
    >>> from skysearcher.skysearch_lib import *
    >>> r = 10
    >>> annuli, annuli_step = get_annuli(r)
    >>> annuli
    [(-3.1730085801256909, -3.1642795331377247),
     (-3.1642795331377247, -3.1555504861497585),
     (-3.1555504861497585, -3.1468214391617924),
     ...
     (3.1468214391617924, 3.1555504861497585),
     (3.1555504861497585, 3.1642795331377247),
     (3.1642795331377247, 3.1730085801256909)]
    >>> annuli_step
    0.0086340369527229677

    """

    _num = ANNULUS_PHI_STEP + int(_r * 1.5)
    # print("num:", _num, "sections |", "radius:", _r, "Kpc")
    arr_, step_ = np.linspace(
        start=-np.pi * 1.01,
        stop=np.pi * 1.01,
        num=_num,
        endpoint=False,
        retstep=True,
        dtype=np.float64)

    """
    # Alternative method.
    global _ANNULUS_PHI_STEP
    if _ANNULUS_PHI_STEP is None:
        _nada, _ANNULUS_PHI_STEP = np.linspace(
            start=-np.pi,
            stop=np.pi,
            num=ANNULUS_PHI_STEP,
            endpoint=False,
            retstep=True,
            dtype=np.float64)
    t1 = (_r + R_STEP)**2 - _r**2
    t2 = (R_START + R_STEP)**2 - R_START**2
    _num = int(((2.0 * np.pi) / _ANNULUS_PHI_STEP) * (t1 / t2))
    arr_ = np.linspace(
        start=-np.pi,
        stop=np.pi,
        num=_num,
        endpoint=False,
        dtype=np.float64)
    print(_num, "sections at radius:", round(_r, 2), "Kpc")
    return list(zip(arr_, arr_ + _ANNULUS_PHI_STEP)), _ANNULUS_PHI_STEP
    """
    return list(zip(arr_, arr_ + step_)), step_


@jit(cache=True)
def get_idx(_grid, _d0, _d1, _ridx):
    """
    The grid index for grid spaces within this segment (:term:`idx`).

    Parameters
    ----------
    _grid : numpy.ndarray
        The data grid for this :term:`halo` (:term:`grid`).
    _d0 : float
        Starting point.
    _d1 : float
        Ending point.
    _ridx : numpy.ndarray
        An array of indices's representing the corresponding array segments.

    Returns
    -------
    numpy.ndarray

    Examples
    --------
    >>> from skysearcher.skysearch_lib import *
    >>> grids = grid_list()
    >>> grid_info = grids[0]
    >>> halo, grid_fh = grid_info
    >>> grid = load_grid(grid_fh)
    >>> grid
    array([[[   0., 0., 0., ..., 0., 423.55697632, -2.3561945 ],
            [   0., 0., 0., ..., 0., 422.85043335, -2.35786676],
            [   0., 0., 0., ..., 0., 422.14511108, -2.35954452],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 2.35954452 ],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35786676 ],
            [   0., 0., 0., ..., 0., 423.55697632, 2.3561945  ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -2.35452223],
            [   0., 0., 0., ..., 0., 422.14276123, -2.3561945 ],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35787249],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 2.35787249 ],
            [   0., 0., 0., ..., 0., 422.14276123, 2.3561945  ],
            [   0., 0., 0., ..., 0., 422.85043335, 2.35452223 ]
            ],
           [[   0., 0., 0., ..., 0., 422.14511108, -2.35284448],
            [   0., 0., 0., ..., 0., 421.43624878, -2.35451674],
            [   0., 0., 0., ..., 0., 420.72854614, -2.3561945 ],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 2.3561945  ],
            [   0., 0., 0., ..., 0., 421.43624878, 2.35451674 ],
            [   0., 0., 0., ..., 0., 422.14511108, 2.35284448 ]],
            ...,
           [[   0., 0., 0., ..., 0., 422.14511108, -0.7887482 ],
            [   0., 0., 0., ..., 0., 421.43624878, -0.787076  ],
            [   0., 0., 0., ..., 0., 420.72854614, -0.78539819],
            ...,
            [   0., 0., 0., ..., 0., 420.72854614, 0.78539819 ],
            [   0., 0., 0., ..., 0., 421.43624878, 0.787076   ],
            [   0., 0., 0., ..., 0., 422.14511108, 0.7887482  ]
            ],
           [[   0., 0., 0., ..., 0., 422.85043335, -0.78707045],
            [   0., 0., 0., ..., 0., 422.14276123, -0.78539819],
            [   0., 0., 0., ..., 0., 421.43624878, -0.78372031],
            ...,
            [   0., 0., 0., ..., 0., 421.43624878, 0.78372031 ],
            [   0., 0., 0., ..., 0., 422.14276123, 0.78539819 ],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78707045 ]
            ],
           [[   0., 0., 0., ..., 0., 423.55697632, -0.78539819],
            [   0., 0., 0., ..., 0., 422.85043335, -0.78372592],
            [   0., 0., 0., ..., 0., 422.14511108, -0.78204811],
            ...,
            [   0., 0., 0., ..., 0., 422.14511108, 0.78204811 ],
            [   0., 0., 0., ..., 0., 422.85043335, 0.78372592 ],
            [   0., 0., 0., ..., 0., 423.55697632, 0.78539819 ]
           ]], dtype=float32)
    >>> grid.shape
    (600, 600, 7)
    >>> _radii = radii()
    >>> r, r_start, r_stop = _radii[0]
    >>> r
    5
    >>> r_start
    4.75
    >>> r_stop
    5.25
    >>> mu, r_idx = mu_idx(grid, r_start, r_stop)
    >>> mu
    6650.8335
    >>> r_idx
    array([[False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           ...,
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False]], dtype=bool)
    >>> annuli, annuli_step = get_annuli(r)
    >>> annuli
    [(-3.1730085801256909, -3.1642795331377247),
     (-3.1642795331377247, -3.1555504861497585),
     (-3.1555504861497585, -3.1468214391617924),
     ...,
     (3.1468214391617924, 3.1555504861497585),
     (3.1555504861497585, 3.1642795331377247),
     (3.1642795331377247, 3.1730085801256909)]
    >>> annuli_step
    0.0087290469879661367
    >>> _deg0, _deg1 = annuli[0]
    >>> _deg0
    -3.1730085801256909
    >>> _deg1
    -3.1642795331377247
    >>> mu, r_idx2 = mu_idx2(grid, r_idx, _deg0, _deg1)
    >>> mu
    1.0
    >>> r_idx2
    array([[False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           ...,
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False],
           [False, False, False, ..., False, False, False]], dtype=bool)
    >>> idx = get_idx(grid, _deg0, _deg1, r_idx2)
    (array([], dtype=int64), array([], dtype=int64))
    >>> type(idx)
    tuple
    >>> type(idx[0])
    numpy.ndarray

    """
    seg_idx_ = np.nonzero(
        np.logical_and(
            np.logical_and(
                _grid[:, :, 6] >= _d0,
                _grid[:, :, 6] < _d1),
            _ridx))
    return seg_idx_


@jit(cache=True)
def get_xbox(_grid, _idx, _mu):
    """
    The array of grid spaces from this segment's contrast density value (:term:`xbox`).

    Parameters
    ----------
    _grid : numpy.ndarray
        The data grid for this halo (:term:`grid`)
    _idx : numpy.ndarray
    _mu : float
        Mean number of stars per array element

    Returns
    -------
    numpy.ndarray

    Example
    -------
    >>> n_boxes_in_seg = len(idx[0])
    >>> xbox = get_xbox(grid, idx, mu)
    >>> n_stars_here = grid[:, :, 1][idx].sum()

    """
    if _mu == 0.0:
        _mu = 1.0
    return (_grid[:, :, 1][_idx] / _mu) - 1.0


def new_sat_stars(_id_lst):
    """
    Start a dictionary for counting stars per :term:`satid` (:term:`sat_stars`).

    Parameters
    ----------
    _id_lst : list
        List of :term:`satid`'s' for this :term:`halo`

    Returns
    -------
    dict
        Dictionary with :term:`satid`'s for keys and corresponding number of stars for each :term:`satid` counted.

    Examples
    --------
    >>> sat_stars = new_sat_stars(satid_list)
    >>> sat_stars
    {0: 0,
     1: 0,
     2: 0,
     3: 0,
     ...
     110: 0,
     111: 0,
     112: 0}

    """
    fresh_dict = {}
    for _id in _id_lst:
        fresh_dict[_id] = 0
    return fresh_dict

def getTerminalSize():
    import platform
    current_os = platform.system()
    tuple_xy=None
    if current_os == 'Windows':
        tuple_xy = _getTerminalSize_windows()
        if tuple_xy is None:
            tuple_xy = _getTerminalSize_tput()
          # needed for window's python in cygwin's xterm!
    if current_os == 'Linux' or current_os == 'Darwin' or  current_os.startswith('CYGWIN'):
        tuple_xy = _getTerminalSize_linux()
    if tuple_xy is None:
        print("default")
        tuple_xy = (80, 25)      # default value
    return tuple_xy

def _getTerminalSize_windows():
    res=None
    try:
        from ctypes import windll, create_string_buffer

        # stdin handle is -10
        # stdout handle is -11
        # stderr handle is -12

        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
    except:
        return None
    if res:
        import struct
        (bufx, bufy, curx, cury, wattr,
         left, top, right, bottom, maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        sizex = right - left + 1
        sizey = bottom - top + 1
        return sizex, sizey
    else:
        return None

def _getTerminalSize_tput():
    # get terminal width
    # src: http://stackoverflow.com/questions/263890/how-do-i-find-the-width-height-of-a-terminal-window
    try:
        import subprocess
        proc=subprocess.Popen(["tput", "cols"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output=proc.communicate(input=None)
        cols=int(output[0])
        proc=subprocess.Popen(["tput", "lines"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output=proc.communicate(input=None)
        rows=int(output[0])
        return (cols,rows)
    except:
        return None


def _getTerminalSize_linux():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,'1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (env['LINES'], env['COLUMNS'])
        except:
            return None
    return int(cr[1]), int(cr[0])

@jit(cache=True)
def count_strs(_dict, _region, _table):
    """
    Count stars per :term:`satid`.

    Parameters
    ----------
    _dict : dict
        Dictionary for stars per :term:`satid` (:term:`sat_stars`).
    _region : list
        Put all region info into a list (:term:`r_info`).
    _table : astropy.table.Table
        :term:`Local_satid_table`: This is so we don't need to index the whole table every loop of the following for loop (:term:`local_satid_table`).

    Returns
    -------
    dict
        An updated :term:`sat_stars` dictionary with :term:`satid`'s for keys and corresponding number of stars for each satid counted.

    Example
    -------
    >>> sat_stars = new_sat_stars(satid_list)
    >>> local_satid_table = satid_table[np.logical_and(
    ...:   satid_table["Rads"] >= r_start,
    ...:   satid_table["Rads"] < r_stop)]
    >>> r_info = [r_start, r_stop, _deg0, _deg1]
    >>> sat_stars = count_strs(sat_stars, r_info, local_satid_table)

    """
    _r_start, _r_stop, _deg0, _deg1 = _region
    _lims = np.nonzero(np.logical_and(
        _table["Phi"] >= _deg0, _table["Phi"] < _deg1))
    sats = _table["satids"][_lims].tolist()
    for _satid in list(_dict.keys()):
        _dict[_satid] += sats.count(_satid)
    _table.remove_rows(_lims)
    return _dict


def dom_satid(input_dict, _dict, rtn=False):
    """
    The dominate satellite number (:term:`domsat_id`).  Domsats's % of all
    stars (:term:`domsat_purity`).

    Parameters
    ----------
    _dict : dict
        A dictionary for counting stars per :term:`satid` (:term:`sat_stars`).

    Returns
    -------
    tuple
        (_domsat_id, _domsat_per, _standout, n_sats)
        :term:`_domsat_id` : Dominate satellite id.
        :term:`_domsat_per` : Percentage that the _domsat_id represents from all satid's.
        :term:`_standout` : _domsat_per / 0.01
        :term:`n_sats` : Number of staid's present.

    Example
    -------
    >>> domsat_id, domsat_purity, standout, nsats = dom_satid(sat_stars)

    """
    best_sat = (1, 0)  # (nstars, id)
    second_best = (1, 0)
    n_sats = 0

    for _id in list(input_dict.keys()):
        if input_dict[_id]:
            n_sats += 1
        if input_dict[_id] > best_sat[0]:
            second_best = best_sat
            best_sat = (input_dict[_id], _id)
    _dict["input_sats"] = input_dict.copy()
    _dict["best_sat"] = best_sat
    _dict["second_best"] = second_best
    _dict["n_sats_start"] = n_sats

    mean = np.mean([i for i in list(input_dict.values()) if i >= 1.0])
    _dict["mean"] = mean

    for key in list(input_dict.keys()):
        if input_dict[key] < mean:
            del input_dict[key]

    _dict["n_sats_end"] = len(input_dict)
    _dict["output_sats"] = input_dict.copy()

    """
    if len(list(_dict.values())) > 20:
        dom_satid(_dict, rtn=True)
    if rtn:
        return
    """

    total_strs = np.sum(list(input_dict.values()))
    _domsat_id = int(best_sat[1])
    _domsat_per = best_sat[0] / float(total_strs)
    _dict["dsat-nstars"] = total_strs
    _dict["domsat_id"] = _domsat_id
    _dict["domsat_purity"] = _domsat_per

    if second_best[0] > 0:
        _standout = _domsat_per / (second_best[0] / float(total_strs))
    else:
        _standout = _domsat_per / 0.01
    _dict["standout"] = _standout

    return _dict

def br(_str, _sclr=CYAN, _bclr=GREEN):

    return ct("[ ", _bclr) + ct(_str, _sclr) + ct(" ]", _bclr)

def prnk(_rnk, _mode="normal"):

    if _mode == "normal":

        p0 = ""
        p1 = ct("[ ", GREEN)
        p2 = ct("MPI_RANK ", CYAN)
        p3 = ct(str(_rnk), MAGENTA)
        p4 = ct(" ]", GREEN)

        return p0 + p1 + p2 + p3 + p4

def report(_type, _info):
    """
    _dict = sorted(_dict, key=_dict.__getitem__, reverse=True)
    if _type == "starting halo":
    STDOUT.write("rank " + str(MPI_RANK) +
    " [ LOADED ] " + halo + "\\n")
    STDOUT.flush()
    if _type == "end annulus":
    line = ("rank " + str(MPI_RANK) + " : halo" + halo[-2:] + " - " +
    str(round(r, 1)) + " Kpc : --> " +
    str(round((time() - a_tic), 1)) + " secs")
    STDOUT.write(line + "\\n")
    STDOUT.flush
    if _type == "end halo":
    if _type == "exit":
    msg = ("rank " + str(MPI_RANK) + " [ FINISHED ] [ " +
    str(round((time() - tic) / 60.0, 1)) + " minutes ]\\n")
    STDOUT.write(msg)
    STDOUT.flush()
    --> TIME
    """

    boxsize, boxheight = getTerminalSize()
    begin_line = ct("| ", WHITE)
    sep_1 = " : "
    sep_2 = " - "
    end_line = ct(" |", WHITE)
    base_len = 10

    def headder(_spc=False, _styl="=", _clr=WHITE):
        if _spc:
            _h1 = "\n\n"
        else:
            _h1 = ""
        _h2 = ct(_styl * boxsize, _clr)
        print(_h1 + _h2)
        STDOUT.flush()

    def divider(_clr=WHITE):
        _size = boxsize - 4
        _1 = ct("| ", _clr)
        _2 = ct("-" * _size, _clr)
        _3 = ct(" |", _clr)
        print(_1 + _2 + _3)
        STDOUT.flush()

    def emptyline(_clr=WHITE):
        _size = boxsize - 4
        _1 = ct("| ", _clr)
        _2 = " " * _size
        _3 = ct(" |", _clr)
        print(_1 + _2 + _3)
        STDOUT.flush()

    # Convert all values to strings.
    for key in list(_info.keys()):
        _dtyp = type(_info[key])
        if _dtyp == str:
            continue
        elif _dtyp in [float, np.float64, np.float32]:
            _info[key] = str(round(_info[key], 2))
        elif _dtyp in [int, dict, tuple, np.int64]:
            _info[key] = str(_info[key])
        else:
            print(ct(_dtyp, RED))
            print(ct(_info[key], RED))
            _info[key] = str(_info[key])

    if _type == "new region":

        # Keys to ignore.
        printed = [
            "input_sats",
            "output_sats",
            "halo_num",
            "MPI_RANK",
            "halo"
        ]

        # Value descriptions.
        defs = {
            "best_sat": "(nstars, satid)",
            "second_best": "(nstars, satid)",
            "n_sats_start": "number of satids before removing < mean",
            "mean": "The mean number of stars per nonzero satid",
            "n_sats_end": "number of satids after removing < mean",
            "dsat-nstars": "total number of stars for all satids",
            "domsat_id": "satid of dominate sat",
            "domsat_purity" :"% (nstars / dsat-nstars)",
            "standout": "% ( domsat_purity / (second_best n_stars / dsat-nstars) )",
            "r_start": "starting radius Kpc",
            "r_stop": "ending radius Kpc",
            "annuli_step": "angular step in Radians",
            "log(mu_2)": "mean number of stars in annulus sub section",
            "max xbox": "max of the max xbox values",
            "max log(nstr)": "max of the log(n_stars)",
            "n_stars": "total stars in region",
            "n_segments": "number of search segments in region",
            "n_seg_increase": "",
            "n_seg_decrease": "",
            "starting deg": "Radians",
            "angular extent": "Radians",
            "halo": "",
            "mass": "M_sun",
            "atime": "GigaYr",
            "jcirc": ""
        }

        # Find column spacing.
        key_lnth = 0
        value_lnth = 0
        for key, value in _info.items():
            if key in printed:
                continue
            if len(key) > key_lnth:
                key_lnth = len(key)

            if len(value) > value_lnth:
                value_lnth = len(value)

        base_len += value_lnth + key_lnth


        # Top of box.
        os.system("clear")
        STDOUT.flush()
        headder()
        divider()

        # Header message.
        _spc = " - "
        _now = datetime.datetime.now()
        _vals = [GREEN, BLUE, RED, YELLOW, GREEN]
        _h1 = _info["halo"].title() + _spc
        _h2 = _type.title() + _spc
        _h3 = _now.strftime("%I:%M:%S %p - %A")
        _h4 = ["[ ", _h1, _h2, _h3, " ]"]
        _size = boxsize - 13 - len(_info["MPI_RANK"])
        _1 = begin_line
        _2 = render_and_pad(_size, zip(_vals, _h4))
        _3 = prnk(_info["MPI_RANK"])
        _4 = end_line
        print(_1 + _2 + _3 + _4)
        STDOUT.flush()

        # Top section of box.
        # Divider.
        divider()

        # Main contents.
        values = [WHITE, GREEN, YELLOW, CYAN, YELLOW, BLUE]
        for key in list(_info.keys()):
            if key in printed:
                continue

            key_space = key_lnth - len(key)
            value_space = value_lnth - len(_info[key])
            _l1 = begin_line
            _l2 = key + " " * key_space
            _l3 = sep_1
            _l4 = _info[key] + " " * value_space
            _l5 = sep_2
            _l6 = defs[key]
            _ln = [_l1, _l2, _l3, _l4, _l5, _l6]
            _size = boxsize  + 14
            _line = zip(values, _ln)
            _1 = render_and_pad(_size, _line)
            _2 = end_line
            print(_1 + _2)
            STDOUT.flush()

        # Midsection of box.
        # Divider.
        divider()

        # Break up satid counts.
        strdict = _info["output_sats"].split(",")
        _minsz = 5
        for x in strdict:
            if len(x) > _minsz:
                _minsz = len(x)
        csize = boxsize // _minsz
        chunks = [strdict[x:x + csize + 4] for x in range(0, len(strdict), csize + 4)]

        # Main content.
        for ln in chunks:
            _line = ",".join(ln)
            _size = boxsize - 4
            _1 = begin_line
            _2 = render_and_pad(_size, [(MAGENTA, _line)])
            _3 = end_line
            print(_1 + _2 + _3)
            STDOUT.flush()

        # Midsection of box.
        # Divider.
        divider()

        # Break up satid counts.
        strdict = _info["input_sats"].split(",")
        chunks = [strdict[x:x + csize] for x in range(0, len(strdict), csize)]

        # Main content.
        for ln in chunks:
            _line = ",".join(ln)
            _size = boxsize - 4
            _1 = begin_line
            _2 = render_and_pad(_size, [(BLUE, _line)])
            _3 = end_line
            print(_1 + _2 + _3)
            STDOUT.flush()

        # Bottom of box.
        headder()

    if _type == "save":
        _color = GREEN
        _vals = [RED, MAGENTA, RED]
        # headder(_styl="-", _clr=_color)
        _l1 = "\t [ "
        _l2 = "SAVING TABLE"
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = "" # ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = "" # ct(" |", _color)
        print(_1 + _2 + _3 + _4)
        # headder(_styl="-", _clr=_color)
        STDOUT.flush()

    if _type == "starting halo":
        _color = BLACK
        _vals = [YELLOW, BLUE, YELLOW]

        _l1 = "\t [ "
        _l2 = "STARTING " + _info["halo"].upper()
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = "" # ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = "" # ct(" |", _color)
        # headder(_clr=_color)
        # emptyline(_clr=_color)
        print(_1 + _2 + _3 + _4)
        # emptyline(_clr=_color)
        # headder(_clr=_color)
        STDOUT.flush()

    if _type == "end annulus":
        _color = BLUE
        _vals = [GREEN, CYAN, GREEN]

        _l1 = "\t [ "
        _l2 = "FINISHED ANNULUS : \t" + _info["annulus time"] + " seconds"
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = "" # ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = "" # ct(" |", _color)
        # headder(_styl="-", _clr=_color)
        # emptyline(_clr=_color)
        print(_1 + _2 + _3 + _4)
        # emptyline(_clr=_color)
        # headder(_styl="-", _clr=_color)
        STDOUT.flush()

    if _type == "end halo":
        _color = MAGENTA
        _vals = [MAGENTA, YELLOW, MAGENTA]
        _l1 = "\t [ "
        _l2 = "FINISHED " + _info["halo"].upper() + " : \t" + _info["halo time"] + " minutes"
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = "" # ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = "" # ct(" |", _color)
        # headder(_styl="-", _clr=_color)
        # emptyline(_clr=_color)
        print(_1 + _2 + _3 + _4)
        # emptyline(_clr=_color)
        # headder(_styl="-", _clr=_color)
        STDOUT.flush()

    if _type == "exit":
        _color = BLACK
        _vals = [RED, RED, RED]
        _l1 = "\t [ "
        _l2 = "FINISHED PROGRAM: \t" + _info["program time"] + " minutes"
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = "" # ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = "" # ct("") |", _color)
        # headder(_clr=_color)
        # emptyline(_clr=_color)
        print(_1 + _2 + _3 + _4)
        # emptyline(_clr=_color)
        # headder(_clr=_color)
        STDOUT.flush()

    if _type == "error":

        # Find column spacing.
        key_lnth = 0
        value_lnth = 0
        for key, value in _info.items():
            if len(key) > key_lnth:
                key_lnth = len(key)

            if len(value) > value_lnth:
                value_lnth = len(value)

        _color = BLACK
        _vals = [YELLOW, RED, YELLOW]
        headder(_styl="-", _clr=_color)
        _l1 = " [ "
        _l2 = "ERROR"
        _l3 = " ]"
        _ln = [_l1, _l2, _l3]
        _size = boxsize - 15 - len(_info["MPI_RANK"])
        _line = zip(_vals, _ln)
        _1 = ct("| ", _color)
        _2 = prnk(_info["MPI_RANK"])
        _3 = render_and_pad(_size, _line)
        _4 = ct(" |", _color)
        print(_1 + _2 + _3 + _4)

        # Main contents.
        values = [_color, RED, YELLOW, RED]
        for key in list(_info.keys()):

            key_space = key_lnth - len(key)
            value_space = value_lnth - len(_info[key])
            _l1 = ct("| ", _color)
            _l2 = key + " " * key_space
            _l3 = sep_1
            _l4 = _info[key] + " " * value_space
            _ln = [_l1, _l2, _l3, _l4]
            _size = boxsize  + 14
            _line = zip(values, _ln)

            _1 = render_and_pad(_size, _line)
            _2 = ct(" |", _color)
            print(_1 + _2)
            STDOUT.flush()

        headder(_styl="-", _clr=_color)
        STDOUT.flush()

    STDOUT.flush()


if __name__ == "__main__":

    now = datetime.datetime.now()
    print(br(now.strftime("%I:%M:%S %p - %A %B %d %Y")))