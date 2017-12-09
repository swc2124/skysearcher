# -*- coding: utf-8 -*-
# @Author: sol courtney
# @Date:   2017-08-23 13:47:59
# @Last Modified by:   swc21
# @Last Modified time: 2017-08-29 23:30:02

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from time import sleep
import configparser as ConfigParser
import os
import pickle
import sys

from random import shuffle

import numpy as np

from numba import jit

from astropy.table import Table
from astropy.table import vstack
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt


# =============================================================================
# =============================================================================

_table_columns = [
    ('radius', 'i'),
    ('r0', 'f'),
    ('r1', 'f'),
    ('annuli_step', 'f'),

    ('halo', 'i'),

    ('xbox_cut', 'f'),
    ('_xbox_fixed', 'f'),
    ('_xbox_mu', 'f'),
    ('Log10(mu)', 'f'),

    ('xbox_min', 'f'),
    ('xbox_mean', 'f'),
    ('xbox_max', 'f'),
    ('Log10(n_stars_max)', 'f'),
    ('Log10(n_stars)', 'f'),

    ('domsat_purity', 'f'),
    ('domsat_id', 'i'),
    ('domsat_mass', 'f'),
    ('domsat_atime', 'f'),

    ('deg0', 'f'),
    ('deg1', 'f'),
    ('extent', 'f'),
    ('n_segments', 'i'),
    ('n_boxes', 'i'),

    ('rank', 'i')]

_badwords = ['brown', 'light', 'grey', 'gray', 'white', 'ivory',
             'beige', 'yellow', 'lime', 'gold', 'salmon', 'khaki',
             'pale', 'olive', 'pink', 'rose', 'silver', 'peach',
             'sea', 'linen', 'coral', 'fuchsia', 'sky', 'tan',
             'snow', 'lemon', 'old', 'mint', 'cream', 'bisque',
             'papaya']


def _make_colors(_clrs=list(mcolors.cnames.keys()), _bdwrds=_badwords):
    for _i, _name in enumerate(_clrs):
        for _bad_word in _bdwrds:
            if _bad_word in _name:
                del _clrs[_i]
                break
    shuffle(_clrs)
    shuffle(_clrs)
    return _clrs

_colors = _make_colors()
_colors = _make_colors(_colors)
_colors = _make_colors(_colors)
_colors = _make_colors(_colors)
_colors = _colors * 150

# print messages
_msg_exit_wrongdir = (
    'Please change to the correct directory before starting'
    'cd skysearcher/skysearcher'
    'Thank you :)')

# select new configuration file msg prompt
_msg_rcfile_select = (
    'select from the following existing configuration files:')
_msg_rcfile_select_line_small = '-' * len(_msg_rcfile_select)
_msg_rcfile_select_line_big = '=' * len(_msg_rcfile_select)
_msg_rcfile_select_rcfile_choice_error = (
    'We are making you a new rc file because you choose wrong.')

# =============================================================================
# =============================================================================


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

# =============================================================================
# =============================================================================

# <configuration file values>
# Make a new configuration file parser object (config).
_config = ConfigParser.ConfigParser()

# Read in rc.file and unpack all values.
_config.read(sortout_rcfile())

# PATH
DATA_DIR = _config.get('PATH', 'data_dir')
GRID_DIR = _config.get('PATH', 'GRID_DIR')
TABLE_DIR = _config.get('PATH', 'table_dir')
MPI_TABLE_DIR =  _config.get('PATH', 'mpi_table_dir')
PLOT_DIR = _config.get('PATH', 'plot_dir')
T_FRMT = _config.get('PATH', 'table_format')
H5_PTH = _config.get('PATH', 'table_hdf5_path')
TBL_EXT = _config.get('PATH', 'table_ext')

# Data
D_MPC = _config.getfloat('Data', 'd_mpc')

# Search Extent
R_START = _config.getint('Search Extent', 'r_start')
R_STOP = _config.getint('Search Extent', 'r_stop')
R_STEP = _config.getint('Search Extent', 'r_step')
R_SCALE = _config.getfloat('Search Extent', 'annulus_scale')
DEG_0 = _config.getfloat('Search Extent', 'a0')
DEG_1 = _config.getfloat('Search Extent', 'a1')
ANNULUS_PHI_STEP = _config.getint('Search Extent', 'annulus_phi_step')

# EXPERIMENTAL PARAMS
EXP_ANGULAR_EXTENT = .2    # percentage of the anulus to evaluate
                            # centered on the current segment.
MIN_LOG_NSTARS = _config.getfloat('Accept Reject', 'min_log_nstars')
MIN_SEG_SIZE = _config.getint('Accept Reject', 'min_seg_size')
MAX_SEG_SIZE = _config.getint('Accept Reject', 'max_seg_size')
MIN_BOXES = _config.getint('Accept Reject', 'min_boxes')
MAX_BOXES = _config.getint('Accept Reject', 'max_boxes')
XBOX_MIN_MU = _config.getfloat('Accept Reject', 'xbox_min_mu')
XBOX_MIN_FIXED = _config.getfloat('Accept Reject', 'xbox_min_fixed')
XBOX_MIN = _config.getfloat('Accept Reject', 'xbox_min')

DIR_LIST = [DATA_DIR, GRID_DIR, TABLE_DIR, PLOT_DIR]

# =============================================================================
# =============================================================================
# Functions bellow require definitions above.

def sortout_directories(directory_list=DIR_LIST):
    '''
    Helper function to quickly check that all directories are in place.

    Keyword Arguments:
        directory_list {list} -- list of directories to check
        (default: {dir_list})

    Returns:
        bool -- False if there were failures
        bool -- True if no failures.
    '''
    # Make an empty list to store failed directory names in.
    failed_list = []
    # Iterate through the given list
    for directory in directory_list:
        # If the directory doesn't exist:
        if not os.path.isdir(directroy):
            # Use the try clause to be safe.
            try:
                # Try to make a directory.
                os.mkdir(directory)
            except:
                # If Mkdir() fails, store name
                failed_list.append(directory)
    # If there were any fails.
    if len(failed_list):
                # Print out each directory name.
        for name in failed_list:
            print(name, 'failed to make directory')
        # return False
        return False
    # No failures, return True
    return True

# =============================================================================
# =============================================================================

def message(_table, _halo, _radius, _xbmin, _mu, _r0, _r1):
    '''Simple helper function for clean stdout.

    Used by skysearcher for consistent output.

    Arguments:
        _table {astropy.table.Table} -- astropy table object
        _halo {str} -- name of the halo i.e. "halo02"
        _radius {int} -- radial value in Kpc
        _xbmin {float} -- the minimum xbox value
        _mu {float} -- the average xbox for the radius
        _r0 {float} -- inner radial value Kpc
        _r1 {float} -- outer radial value Kpc
    '''
    if os.name == 'nt':
        clr_cmd = 'cls'
    else:
        clr_cmd = 'clear'
    os.system(clr_cmd)
    r = str(_radius)
    msg0 = _halo + '    -   radius: ' + r + ' Kpc   -'
    msg00 = 'r0:' + str(_r0) + '  r1:' + str(_r1) + '   -'
    msg1 = 'xbox_min: ' + str(round(_xbmin, 3)) + '     -'
    msg2 = 'mu: ' + str(round(_mu, 3)) + ' strs/kpc^2   -'
    print(msg0 + msg1 + msg2)
    if len(_table):
        _table.pprint(max_lines=500, max_width=200, align='^')

@jit(cache=True)
def xbox(_grid, _idx, _mu):
    return (_grid[:, :, 0][_idx] / _mu) - 1.0

#######################################################################################
# -- UNDER CONSTRUCTION ------------------------------------------------------------- #
#######################################################################################

@jit(cache=True)
def EXP_mu_idx(_grid, _d0, _d1, r0, r1):
    if _d0 < -np.pi or _d1 > np.pi:
        #print(_d0, _d1)
        grid_idx_ = np.logical_and(
                        np.logical_and(
                            _grid[:, :, 4] >= r0,
                            _grid[:, :, 4] < r1),
                        np.logical_and(
                            _grid[:, :, 5].T >= (-_d0 + (1.5 * np.pi)),
                            _grid[:, :, 5].T < (-_d1 + (1.5 * np.pi))))       

    else:
        grid_idx_ = np.logical_and(
                        np.logical_and(
                            _grid[:, :, 4] >= r0,
                            _grid[:, :, 4] < r1),
                        np.logical_and(
                            _grid[:, :, 5] >= _d0,
                            _grid[:, :, 5] < _d1))
    return _grid[:, :, 0][np.nonzero(grid_idx_)].mean(), grid_idx_

#######################################################################################
# ----------------------------------------------------------------------------------- #
#######################################################################################

@jit(cache=True)
def kpc_to_arcmin(d_mpc=D_MPC):
    d = 1e3  # Kpc
    D = d_mpc * 1e6
    arcmin_mod = 3437.746  # (60.0 * 180.0) / np.piS
    return np.square(arcmin_mod * (d / D))

@jit(cache=True)
def get_idx(_grid, _d0, _d1, _ridx):
    seg_idx_ = np.nonzero(
        np.logical_and(
            np.logical_and(
                _grid[:, :, 5] >= _d0,
                _grid[:, :, 5] < _d1),
            _ridx))
    return seg_idx_

def grid_list():
    return [(fh.split('_')[0], os.path.join(GRID_DIR, fh))
            for fh in os.listdir(GRID_DIR)
            if fh.endswith('npy')]

def load_grid(_grid_fh):
    return fix_rslice(np.load(_grid_fh))

@jit(cache=True)
def mu_idx(_grid, r0, r1):
    grid_idx_ = np.logical_and(
        _grid[:, :, 4] >= r0,
        _grid[:, :, 4] < r1)
    return _grid[:, :, 0][np.nonzero(grid_idx_)].mean(), grid_idx_

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
        r_tbl_.meta['deg_0'] = DEG_0
        r_tbl_.meta['deg_1'] = DEG_1
        r_tbl_.meta['min_seg_size'] = MIN_SEG_SIZE
    return r_tbl_

def radii():
    radii_list_ = []
    for r in list(range(R_START, R_STOP, R_STEP)):
        annulus_scale_ = R_SCALE * r
        if r >= 100:
            annulus_scale_ *= 0.6
            #pass
        radii_list_.append((r, r - annulus_scale_, r + annulus_scale_))
    return radii_list_

def annuli(_r):
    _num = ANNULUS_PHI_STEP + (2 * _r)
    print('num:', _num, 'sections |', 'radius:', _r, 'Kpc')
    arr_, step_ = np.linspace(
                        start=DEG_0,
                        stop=DEG_1,
                        num=_num,
                        endpoint=False,
                        retstep=True,
                        dtype=np.float32)
    return list(zip(arr_, arr_ + step_)), step_

@jit(cache=True)
def fix_rslice(_grid_):
    _grid_ = np.concatenate((_grid_, np.zeros((601, 601, 1))), axis=2)
    center = 301.0
    for i in range(_grid_.shape[0]):
        for q in range(_grid_.shape[1]):
            value = np.sqrt(
                np.square(i - center) + np.square(q - center))
            _grid_[i, q, 4] = value
            _grid_[i, q, 5] = np.arctan2(q - center, i - center)
    return _grid_

def dom_satid(_dict):
    best_sat = (0, 0)  # (nstars, id)
    total_strs = np.sum(list(_dict.values()))
    if not total_strs:
        return False, False
    for _id in _dict.keys():
        if _dict[_id] > best_sat[0]:
            best_sat = (_dict[_id], _id)
    try:
        return int(best_sat[1]), best_sat[0] / (1e0 * total_strs)
    except ZeroDivisionError as e:
        sys.stdout.write(e[0] + '\n')
        return False, False

def satid_setup(_halo, attempts=0):
    table_fh = os.path.join(TABLE_DIR, _halo + TBL_EXT)
    try:
        _tbl = Table.read(table_fh, format=T_FRMT, path=H5_PTH)
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

def count_strs(_dict, _region, _table):
    try:
        R_START, R_STOP, _deg0, _deg1 = _region
        _lims = np.nonzero(
            np.logical_and(
                np.logical_and(
                    _table['Phi'] >= _deg0,
                    _table['Phi'] < _deg1),
                np.logical_and(
                    _table['Rads'] >= R_START,
                    _table['Rads'] < R_STOP)))
        sats = _table['satids'][_lims].tolist()
        for _satid in list(_dict.keys()):
            _dict[_satid] += sats.count(_satid)
        _table.remove_rows(_lims)
        return _dict
    except MemoryError as e:
        print(e)
        return _dict

def new_sat_stars(_id_lst):
    fresh_dict = {}
    for _id in _id_lst:
        fresh_dict[_id] = 0
    return fresh_dict

def update_plot(_ax, _position, _halo, _dom_satid):
    # plt.bar(left, height, width, bottom, hold=None, data=None)
    theta0, r_extent, angular_extent, R_START = _position
    # Plot the regions area.
    _bump = 0.5 * angular_extent
    theta0 += _bump
    _ax.bar(theta0, r_extent,
            angular_extent, R_START,
            color=_colors[_dom_satid],
            alpha=.2)
    # Update/save the figure.
    fig = _ax.get_figure()
    plot_fh = os.path.join(PLOT_DIR, _halo + '_dataplot.png')
    fig.savefig(plot_fh)

def save_plot(_ax, _halo, _out_dir=PLOT_DIR):
    fig = _ax.get_figure()
    plot_fh = os.path.join(_out_dir, _halo + '_dataplot.png')
    fig.savefig(plot_fh)

def plot_full_halo(_halo, _d_cut=10, _out_dir=PLOT_DIR):
    fig = plt.figure(figsize=(20, 20))
    ax1 = fig.add_subplot(111, projection='polar')
    ax1.set_title(_halo)
    table_fh = os.path.join(TABLE_DIR, _halo + TBL_EXT)
    table = Table.read(table_fh, format=T_FRMT, path=H5_PTH)
    ax1.scatter(table['Phi'][::10], table['Rads'][::10],
                s=10,
                alpha=.075,
                marker='.',
                cmap=plt.cm.Paired,
                c=table['Xbox'][::10],
                vmin=0.0,
                vmax=12.0)
    ax1.set_thetagrids(np.linspace(0.0, 360.0, 32)[:-1])
    ax1.set_theta_direction(-1)
    ax1.set_theta_zero_location("N")
    ax1.set_ylim([0, 300])
    plot_fh = os.path.join(_out_dir, _halo + '_fullplot.png')
    fig.savefig(plot_fh)
    plt.close()

def get_data_plot(_halo):
    fig = plt.figure(figsize=(20, 20))
    ax1 = fig.add_subplot(111, projection='polar')
    ax1.set_title(_halo)
    ax1.set_thetagrids(np.linspace(0.0, 360.0, 32)[:-1])
    ax1.set_theta_direction(-1)
    ax1.set_theta_zero_location("N")
    ax1.set_ylim([0, 300])
    return ax1

def log_satstars(_dict, _halo, _r, _dom_satid):
    space = 8
    fh = ('_').join([
        _halo,
        str(int(_r)) + 'kpc',
        'sat' + str(_dom_satid)
        ]) + '.txt'
    satid_log_fh = os.path.join(DATA_DIR, 'satid_logs', fh)
    t_stars = sum(_dict.values())
    if not t_stars:
        return
    with open(satid_log_fh, 'w') as satid_log:
        satid_log.write( 'satid   ' + 'n_stars ' + 'percent\n')
        for key in _dict.keys():
            satid = str(key)
            n_stars = str(_dict[key])
            percent = str(round(1e2 *(_dict[key] / t_stars), 2)) + ' %'
            for i in range(space - len(satid)):
                satid += ' '
            for i in range((space) - len(n_stars)):
                n_stars += ' '
            for i in range(space - len(percent)):
                percent += ' '
            satid_log.write( satid + n_stars + percent + '\n')

def mass_book():
    m_arr = np.load(os.path.join(DATA_DIR, 'satmass_array.npy'))
    _mdict = {}
    for sat in m_arr:
        _mdict[int(sat[0]) - 1] =  round(sat[1], 3)
    return _mdict