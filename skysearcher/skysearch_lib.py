from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from time import sleep
from time import time
import configparser
import os
import pickle
import sys

from mpi4py import MPI
import numpy as np

from numba import jit

from astropy.table import Table
from astropy.table import vstack

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

stdout = os.sys.stdout

def sortout_rcfile():

    #   configuration file settings
    rcfile_ext = os.path.extsep + 'cfg'
    new_rc_fh = 'rc' + rcfile_ext

    # ====================================+====================================
    #                               USER SETTINGS
    # ====================================+====================================
    # <Make sure we're in skysearcher/skysearcher & save path variables>
    # --------------------------------------------------------------------
    # Get a string value for the path to the current directory
    # (curdir_path).
    curdir_path = os.path.abspath(os.path.curdir)

    # Get a string value for
    # 1.) path to parent directory (pardir_path) &
    # 2.) name of current directory (curdir_name).
    pardir_path, curdir_name = os.path.split(curdir_path)

    # Get a string value for parent directory name (pardir_name).
    pardir_name = os.path.split(pardir_path)[1]

    # Make sure that parent & current directory names are the same &
    # if not, then ask user to quit and start again from the proper
    # directory skysearcher/skysearcher.
    if not curdir_name == pardir_name:
        os.sys.exit(_msg_exit_wrongdir)

    # <get configuration file either made or located and then loaded>
    # --------------------------------------------------------------------
    # Determine if there exists a configuration file by making a list
    # of filenames for all the files in the current directory
    # (curdir_files).
    curdir_files = os.listdir(curdir_path)

    # Isolate all files with .cfg extension within curdir_files
    # (rc_files).
    rc_files = [f for f in curdir_files if os.path.splitext(f)[
        1] == rcfile_ext]

    # Ask if there is at least one configuration file present?
    # If there is:
    if len(rc_files):

        # Ask if there is more than one configuration file.
        # If there is:
        if len(rc_files) > 1:
            # Print the top of the prompt
            print(_msg_rcfile_select_line_big)
            print(_msg_rcfile_select)
            print(_msg_rcfile_select_line_small, '\n')
            # Print out available options.
            for i, rc_filename in enumerate(rc_files):
                print('     [', i, ']   ', rc_filename)
            # Ask user for selection.
            selection = int(raw_input('\n\n       Enter selection: '))
            try:
                # Try to use the user input.
                rc_fh = rc_files.pop(selection)
            except IndexError:
                # Finish prompt.
                # Make a new configuration file and clean up.
                import new_cfg
                print(_msg_rcfile_select_rcfile_choice_error)
                rc_fh = new_cfg.new_rc(new_rc_fh)
            finally:
                # Finish prompt and cleanup.
                print(_msg_rcfile_select_line_big)
        # If there is only one:
        else:
            # Take the only name in the list.
            rc_fh = rc_files.pop()
    # If there is not:
    else:
        # Use the configuration file generator function new_cfg.new_rc()
        # and clean up.
        import new_cfg
        rc_fh = new_cfg.new_rc(new_rc_fh)
    return rc_fh

# <configuration file values>
# Make a new configuration file parser object (config).
_config = configparser.ConfigParser()

# Read in rc.file and unpack all values.
_config.read(sortout_rcfile())
D_MPC = _config.getfloat('Data', 'd_mpc')
TABLE_DIR = _config.get('PATH', 'table_dir')
PLOT_DIR = _config.get('PATH', 'plot_dir')
MPI_TABLE_DIR = _config.get('PATH', 'mpi_table_dir')
TABLE_EXT =  _config.get('PATH', 'table_ext')
TABLE_FORMAT =  _config.get('PATH', 'table_format')
TABLE_HDF5_PATH =  _config.get('PATH', 'table_hdf5_path')
GRID_DIR = _config.get('PATH', 'grid_dir')
GRID_EXT =  _config.get('PATH', 'grid_ext')
R_START = _config.getint('Search Extent', 'r_start')
R_STOP = _config.getint('Search Extent', 'r_stop')
R_STEP = _config.getint('Search Extent', 'r_step')
R_SCALE = _config.getfloat('Search Extent', 'annulus_scale')
ANNULUS_PHI_STEP = _config.getint('Search Extent', 'annulus_phi_step')
XBOX_CUT = _config.getfloat('Accept Reject', 'xbox_cut')
N_SKIPS = _config.getint('Accept Reject', 'n_skips')
MIN_LOG_NSTARS = _config.getfloat('Accept Reject', 'min_log_nstars')
MIN_N_SEGMENTS = _config.getint('Accept Reject', 'min_n_segments')
SAVE_INTERVAL = _config.getint('Run Time', 'save_interval')

def pause(rnk=rank):
    sleep(rank * 0.005)

def grid_list():
    return [(fh.split('_')[0], os.path.join(GRID_DIR, fh))
            for fh in os.listdir(GRID_DIR)
            if fh.endswith(GRID_EXT)]

@jit(cache=True)
def kpc_to_arcmin(d_mpc=D_MPC):
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    arcmin_mod = 3437.746  # (60.0 * 180.0) / np.piS
    return np.square(arcmin_mod * (d / D))

_table_columns = [

    # Halo.
    ('halo', 'i'),

    # Annulus location values.
    ('radius', 'i'),
    ('r0', 'f'),
    ('r1', 'f'),
    ('annuli_step', 'f'),
    # ('n_empty_segments', 'i'),

    # Annulus content values.
    # ('log_n_boxes_in_ann', 'i'),
    # ('log_n_stars_in_ann', 'i'),
    ('Log10(mu)', 'f'),

    # Feature content values.
    ('xbox_max', 'f'),
    ('Log10(n_stars_max)', 'f'),
    ('domsat_purity', 'f'),
    ('domsat_id', 'i'),
    ('domsat_mass', 'f'),
    ('domsat_atime', 'f'),

    # Feature location values.
    ('deg0', 'f'),
    # ('deg1', 'f'),
    ('extent', 'f'),
    # ('n_segments', 'i'),
    # ('n_boxes', 'i'),

    # ('rank', 'i')
    ]

def record_table(_names=_table_columns, _meta=True):
    output_colums_ = []
    output_dtyps_ = []
    for name_, dtyp_ in _names:
        output_colums_.append(name_)
        output_dtyps_.append(dtyp_)
    r_tbl_ = Table(names=output_colums_, dtype=output_dtyps_)
    if _meta:
        r_tbl_.meta['r_start'] = R_START
        r_tbl_.meta['r_stop'] = R_STOP
        r_tbl_.meta['r_step'] = R_STEP
        r_tbl_.meta['r_scale'] = R_SCALE
        r_tbl_.meta['fh'] = os.path.join(
                                    MPI_TABLE_DIR,
                                    str(rank) + TABLE_EXT)
    return r_tbl_

def radii():
    radii_list_ = []
    for r in list(range(R_START, R_STOP, R_STEP)):
        annulus_scale_ = R_SCALE * r
        if r >= 75:
            annulus_scale_ *= 0.6
            #pass
        radii_list_.append((r, r - annulus_scale_, r + annulus_scale_))
    return radii_list_

def save_record_table(_table, _rnk=rank):
    _table.write(_table.meta['fh'],
                 format=TABLE_FORMAT,
                 path=TABLE_HDF5_PATH,
                 overwrite=True,
                 append=True)
    stdout.write('rank ' + str(_rnk) + ' [ SAVED ]\n')
    stdout.flush()


@jit(cache=True)
def fix_rslice(_grid):
    grid_ = np.concatenate((_grid.astype(np.float64), np.zeros((601, 601, 1), dtype=np.float64)), axis=2)
    center = 300.0
    for i in range(grid_.shape[0]):
        for q in range(grid_.shape[1]):
            value = np.sqrt(
                np.square(i - center, dtype=np.float64) + np.square(q - center, dtype=np.float64), 
                dtype=np.float64)
            grid_[i, q, 4] = value
            grid_[i, q, 5] = np.arctan2(q - center, i - center)
    return grid_

def load_grid(_grid_fh):

    return fix_rslice(np.load(_grid_fh))

def satid_setup(_halo, attempts=0):
    table_fh = os.path.join(TABLE_DIR, _halo + TABLE_EXT)
    try:
        _tbl = Table.read(table_fh, format=TABLE_FORMAT, path=TABLE_HDF5_PATH)
    except MemoryError as e:
        print(e)
        print('attempt:', attempts, 'failed')
        sleep(5)
        if attempts <= 3:
            satid_setup(_halo=_halo, attempts=(attempts+1))
        else:
            sys.exit(0)
            
    keeps = ['satids', 'Rads', 'Phi']
    _tbl.keep_columns(keeps)
    return _tbl.meta['satids'], _tbl

@jit(cache=True)
def mu_idx(_grid, r0, r1):
    grid_idx_ = np.logical_and(
        _grid[:, :, 4] >= r0,
        _grid[:, :, 4] < r1)
    return _grid[:, :, 0][np.nonzero(grid_idx_)].mean(), grid_idx_

def get_annuli(_r):
    _num = ANNULUS_PHI_STEP + int(_r * 1.5)
    print('num:', _num, 'sections |', 'radius:', _r, 'Kpc')
    arr_, step_ = np.linspace(
                        start=-np.pi * 1.01,
                        stop=np.pi * 1.01,
                        num=_num,
                        endpoint=False,
                        retstep=True,
                        dtype=np.float64)
    return list(zip(arr_, arr_ + step_)), step_

@jit(cache=True)
def get_idx(_grid, _d0, _d1, _ridx):
    seg_idx_ = np.nonzero(
        np.logical_and(
            np.logical_and(
                _grid[:, :, 5] >= _d0,
                _grid[:, :, 5] < _d1),
            _ridx))
    return seg_idx_

@jit(cache=True)
def get_xbox(_grid, _idx, _mu):

    return (_grid[:, :, 0][_idx] / _mu) - 1.0

def new_sat_stars(_id_lst):
    fresh_dict = {}
    for _id in _id_lst:
        fresh_dict[_id] = 0
    return fresh_dict

def count_strs(_dict, _region, _table):
    try:
        _r_start, _r_stop, _deg0, _deg1 = _region
        _lims = np.nonzero(
            np.logical_and(
                np.logical_and(
                    _table['Phi'] >= _deg0,
                    _table['Phi'] < _deg1),
                np.logical_and(
                    _table['Rads'] >= _r_start,
                    _table['Rads'] < _r_stop)))
        sats = _table['satids'][_lims].tolist()
        for _satid in list(_dict.keys()):
            _dict[_satid] += sats.count(_satid)
        _table.remove_rows(_lims)
        return _dict
    except MemoryError as e:
        print(e)
        return _dict

def dom_satid(_dict):
    best_sat = (0, 0)  # (nstars, id)
    total_strs = np.sum(list(_dict.values()))
    for _id in _dict.keys():
        if _dict[_id] > best_sat[0]:
            best_sat = (_dict[_id], _id)
    return int(best_sat[1]), best_sat[0] / (1e0 * total_strs)