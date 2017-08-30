'''[summary]

[description]

Attributes
----------
_table_columns : {list}
    [description]
_badwords : {list}
    [description]
_colors : {[type]}
    [description]
_colors : {[type]}
    [description]
_colors : {[type]}
    [description]
_colors : {[type]}
    [description]
_msg_exit_wrongdir : {tuple}
    [description]
_msg_rcfile_select : {tuple}
    [description]
_msg_rcfile_select_line_small : {str}
    [description]
_msg_rcfile_select_line_big : {str}
    [description]
_msg_rcfile_select_rcfile_choice_error : {tuple}
    [description]
_config : {[type]}
    [description]
_config.read(sortout_rcfile()) : {[type]}
    [description]
data_dir : {[type]}
    [description]
grid_dir : {[type]}
    [description]
table_dir : {[type]}
    [description]
plot_dir : {[type]}
    [description]
t_frmt : {[type]}
    [description]
h5_pth : {[type]}
    [description]
tbl_ext : {[type]}
    [description]
r_start : {[type]}
    [description]
r_stop : {[type]}
    [description]
r_step : {[type]}
    [description]
r_scale : {[type]}
    [description]
deg_0 : {[type]}
    [description]
deg_1 : {[type]}
    [description]
min_seg_size : {[type]}
    [description]
_dir_list : {list}
    [description]
'''
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import ConfigParser
import os
import pickle
import sys

from random import shuffle

import numpy as np

from astropy.table import Table
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt


# =============================================================================
# =============================================================================

_table_columns = [
    ('radius', 'i'),
    ('r0', 'f'),
    ('r1', 'f'),


    ('halo', 'i'),

    ('xbox_cut', 'f'),
    ('_xbox_fixed', 'f'),
    ('_xbox_mu', 'f'),
    ('Log10(mu)', 'f'),

    ('xbox_min', 'f'),
    ('xbox_mean', 'f'),
    ('xbox_max', 'f'),

    ('deg0', 'f'),
    ('deg1', 'f'),
    ('extent', 'f'),
    ('n_segments', 'i'),
    ('arc_len', 'f'),

    ('sat_purity', 'f'),
    ('domsat', 'i'),
    ('domsat_mass', 'f'),

    ('Log10(n_stars)', 'f'),

    ('n_boxes', 'i'),

    ('rank', 'i')]

_badwords = ['brown', 'light', 'grey', 'gray', 'white', 'ivory',
             'beige', 'yellow', 'lime', 'gold', 'salmon', 'khaki',
             'pale', 'olive', 'pink', 'rose', 'silver', 'peach',
             'sea', 'linen', 'coral', 'fuchsia', 'sky', 'tan',
             'snow', 'lemon', 'old', 'mint', 'cream', 'bisque',
             'papaya']


def _make_colors(_clrs=mcolors.cnames.keys(), _bdwrds=_badwords):
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
data_dir = _config.get('PATH', 'data_dir')
grid_dir = _config.get('PATH', 'grid_dir')
table_dir = _config.get('PATH', 'table_dir')
plot_dir = _config.get('PATH', 'plot_dir')
t_frmt = _config.get('PATH', 'table_format')
h5_pth = _config.get('PATH', 'table_hdf5_path')
tbl_ext = _config.get('PATH', 'table_ext')

# Search Extent
r_start = _config.getint('Search Extent', 'r_start')
r_stop = _config.getint('Search Extent', 'r_stop')
r_step = _config.getint('Search Extent', 'r_step')
r_scale = _config.getfloat('Search Extent', 'annulus_scale')
deg_0 = _config.getfloat('Search Extent', 'a0')
deg_1 = _config.getfloat('Search Extent', 'a1')
#deg_step = _config.getfloat('Search Extent', 'annulus_phi_step')
min_seg_size = _config.getint('Search Extent', 'min_seg_size')
xbox_min_mu = _config.getint('Search Extent', 'xbox_min_mu')
xbox_min_fixed = _config.getint('Search Extent', 'xbox_min_fixed')

_dir_list = [data_dir, grid_dir, table_dir, plot_dir]

# =============================================================================
# =============================================================================
# Functions bellow require definitions above.


def sortout_directories(directory_list=_dir_list):
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


def xbox(_grid, _idx, _mu):
    return (_grid[:, :, 0][_idx] / _mu) - 1.0


def get_idx(_grid, _d0, _d1, _ridx):
    seg_idx = np.nonzero(
        np.logical_and(
            np.logical_and(
                _grid[:, :, 5] >= _d0,
                _grid[:, :, 5] < _d1),
            _ridx))
    return seg_idx


def grid_list():
    return [(fh.split('_')[0], os.path.join(grid_dir, fh))
            for fh in os.listdir(grid_dir)
            if fh.endswith('npy')]


def load_grid(_grid_fh):
    return fix_rslice(np.load(_grid_fh))


def mu_idx(_grid, r0, r1):
    _grid_idx = np.logical_and(
        _grid[:, :, 4] >= r0,
        _grid[:, :, 4] < r1)

    return _grid[:, :, 0][np.nonzero(_grid_idx)].mean(), _grid_idx


def record_table(_dict=_table_columns):
    output_colums = []
    output_dtyps = []
    for _name, _dtyp in _table_columns:
        output_colums.append(_name)
        output_dtyps.append(_dtyp)
    r_tbl = Table(names=output_colums, dtype=output_dtyps)
    r_tbl.meta['r_start'] = r_start
    r_tbl.meta['r_stop'] = r_stop
    r_tbl.meta['r_step'] = r_step
    r_tbl.meta['r_scale'] = r_scale
    r_tbl.meta['deg_0'] = deg_0
    r_tbl.meta['deg_1'] = deg_1
    r_tbl.meta['min_seg_size'] = min_seg_size
    return r_tbl


def radii():
    _radii_list = []
    for r in range(r_start, r_stop, r_step):
        annulus_scale = r * r_scale
        #annulus_scale = r_step / 2.0
        if (2.0 * annulus_scale) < r_step:
            annulus_scale = 1.0
        _radii_list.append((r, r - annulus_scale, r + annulus_scale))
    return np.ascontiguousarray(_radii_list)


def annuli(_r):
    annuli_list = []
    _arr, _step = np.linspace(
        deg_0,
        deg_1,
        4 * _r,
        endpoint=False,
        retstep=True)
    for _sgmt in _arr:
        annuli_list.append((_sgmt, _sgmt + _step))
    return np.ascontiguousarray(annuli_list), _step


def fix_rslice(_grid_):
    phi_slice = np.zeros((601, 601, 1))
    _grid_ = np.concatenate((_grid_, phi_slice), axis=2)
    center = 300.0
    for i in range(_grid_.shape[0]):
        for q in range(_grid_.shape[1]):
            value = np.sqrt(
                np.square(i - center) + np.square(q - center))
            _grid_[i, q, 4] = value
            _grid_[i, q, 5] = np.arctan2(q - center, i - center)
    return _grid_


def dom_satid(_dict):
    best_sat = (0, 0)  # (nstars, id)
    total_strs = sum(_dict.values())
    if not total_strs:
        return False, False
    for _id in _dict.keys():
        if _dict[_id] > best_sat[0]:
            best_sat = (_dict[_id], _id)

    try:
        return best_sat[1], best_sat[0] / (1e0 * total_strs)
    except ZeroDivisionError as e:
        sys.stdout.write(e[0] + '\n')
        return False, False


def satid_setup(_halo):
    table_fh = os.path.join(table_dir, _halo + tbl_ext)
    _tbl = Table.read(table_fh, format=t_frmt, path=h5_pth)
    keeps = ['satids', 'Rads', 'Phi']
    _tbl.keep_columns(keeps)
    return _tbl.meta['satids'], _tbl


def count_strs(_dict, _region, _table):
    r_start, r_stop, _deg0, _deg1 = _region
    _lims = np.nonzero(
        np.logical_and(
            np.logical_and(
                _table['Phi'] >= _deg0,
                _table['Phi'] < _deg1),
            np.logical_and(
                _table['Rads'] >= r_start,
                _table['Rads'] < r_stop)))
    sats = _table['satids'][_lims].tolist()
    for _satid in _dict.keys():
        _dict[_satid] += sats.count(_satid)
    _table.remove_rows(_lims)
    return _dict


def new_sat_stars(_id_lst):
    fresh_dict = {}
    for _id in _id_lst:
        fresh_dict[_id] = 0
    return fresh_dict


def update_plot(_ax, _position, _halo, _dom_satid):
    #   plt.bar(left, height, width, bottom, hold=None, data=None)
    theta0, r_extent, angular_extent, r_start = _position
    # Plot the regions area.
    _bump = 0.5 * angular_extent
    theta0 += _bump
    _ax.bar(theta0, r_extent,
            angular_extent, r_start,
            color=_colors[_dom_satid],
            alpha=.2)
    # Update/save the figure.
    fig = _ax.get_figure()
    plot_fh = os.path.join(plot_dir, _halo + '_dataplot.png')
    fig.savefig(plot_fh)


def save_plot(_ax, _halo):
    fig = _ax.get_figure()
    plot_fh = os.path.join(plot_dir, _halo + '_dataplot.png')
    fig.savefig(plot_fh)


def plot_full_halo(_halo, _d_cut=10):
    fig = plt.figure(figsize=(20, 20))
    ax1 = fig.add_subplot(111, projection='polar')
    ax1.set_title(_halo)
    table_fh = os.path.join(table_dir, _halo + tbl_ext)
    table = Table.read(table_fh, format=t_frmt, path=h5_pth)
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
    plot_fh = os.path.join(plot_dir, _halo + '_fullplot.png')
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
    satid_log_fh = os.path.join(data_dir, 'satid_logs', fh)
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
    m_arr = np.load(os.path.join(data_dir, 'satmass_array.npy'))
    _mdict = {}
    for sat in m_arr:
        _mdict[int(sat[0]) - 1] =  round(sat[1], 3)
    return _mdict