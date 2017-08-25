'''
===============================================================================
skysearcher v0.1.0
August 2017
Sol W. Courtney
swc2124@columbia.edu
Columbia University NYC, NY
===============================================================================
python skysearcher/skysearch_lib.py
'''
# does this need to be here?
from __future__ import division, absolute_import, print_function

import os

import ConfigParser

import numpy as np
from astropy.table import Table


# =============================================================================
# =============================================================================

# print messages
msg_exit_wrongdir = (
    'Please change to the correct directory before starting'
    'cd skysearcher/skysearcher'
    'Thank you :)')

# select new configuration file msg prompt
msg_rcfile_select = (
    'select from the following existing configuration files:')
msg_rcfile_select_line_small = '-' * len(msg_rcfile_select)
msg_rcfile_select_line_big = '=' * len(msg_rcfile_select)
msg_rcfile_select_rcfile_choice_error = (
    'We are making you a new rc file because you choose wrong.')

# =============================================================================
# =============================================================================


def sortout_rcfile():

    # configuration file settings
    rcfile_ext = os.path.extsep + 'cfg'
    new_rc_fh = 'rc' + rcfile_ext

    # ====================================+====================================
    #                               USER SETTINGS
    # ====================================+====================================
    # <Make sure we're in skysearcher/skysearcher & save path variables>
    # ------------------------------------------------------------------
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
        os.sys.exit(msg_exit_wrongdir)

    # <get configuration file either made or located and then loaded>
    # ------------------------------------------------------------------
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
            print(msg_rcfile_select_line_big)
            print(msg_rcfile_select)
            print(msg_rcfile_select_line_small, '\n')
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
                print(msg_rcfile_select_rcfile_choice_error)
                rc_fh = new_cfg.new_rc(new_rc_fh)
            finally:
                # Finish prompt and cleanup.
                print(msg_rcfile_select_line_big)

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
config = ConfigParser.ConfigParser()

# Read in rc.file and unpack all values.
config.read(sortout_rcfile())

# Data
data_cut = config.getint('Data', 'data_cut')
halo = config.get('Data', 'halo')
d_mpc = config.getfloat('Data', 'd_mpc')
filter_type = config.get('Data', 'filter_type')

# PATH
curdir_path = config.get('PATH', 'curdir_path')
curdir_name = config.get('PATH', 'curdir_name')
pardir_path = config.get('PATH', 'pardir_path')
pardir_name = config.get('PATH', 'pardir_name')
data_dir = config.get('PATH', 'data_dir')
grid_dir = config.get('PATH', 'grid_dir')
table_dir = config.get('PATH', 'table_dir')
plot_dir = config.get('PATH', 'plot_dir')

t_format = config.get('PATH', 'table_format')
hdf5_pth = config.get('PATH', 'table_hdf5_path')

# Program Units
units_radius = config.get('Program Units', 'units_radius')
units_annulus = config.get('Program Units', 'units_annulus')

# Search Extent
r_start = config.getint('Search Extent', 'r_start')
r_stop = config.getint('Search Extent', 'r_stop')
annulus_scale = config.getfloat('Search Extent', 'annulus_scale')
annulus_phi_step = config.getfloat('Search Extent', 'annulus_phi_step')
xbox_min = config.getfloat('Search Extent', 'xbox_min')

# Accept Reject
min_boxes = config.getint('Accept Reject', 'min_boxes')

# Plots
full_projection = config.getboolean('Plots', 'full_projection')
heatmap = config.getboolean('Plots', 'heatmap')
regions_together = config.getboolean('Plots', 'regions_together')


# Grid file name without full PATH (grid_file_name).
grid_file_name = '_'.join([config.get('PATH', 'fh_prefix'), config.get(
    'PATH', 'grid_file_designator')]) + config.get('PATH', 'grid_ext')


# Grid file handle (full Path) (grid_fh).
grid_fh = os.path.join(config.get('PATH', 'grid_dir'), grid_file_name)

# Table file name without full PATH (table_file_name).
table_file_name = '_'.join([config.get('PATH', 'fh_prefix'), config.get(
    'PATH', 'table_file_designator')]) + config.get('PATH', 'table_ext')


# Table file handle (full Path) (table_fh).
table_fh = os.path.join(config.get('PATH', 'table_dir'), table_file_name)

dir_list = [data_dir, grid_dir, table_dir, plot_dir]


# =============================================================================
# =============================================================================
# Functions bellow require definitions above.


def sortout_directories(directory_list=dir_list):
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


def load_halo(
        g_fh=grid_fh,
        t_fh=table_fh,
        frmt=t_format,
        h5_pth=hdf5_pth):
    '''
    This loads both a grid and table for the halo in the configuration
    file. Just for easy loading, returns both grid and table.

    Keyword Arguments:
        g_fh {str} -- grid file handle (default: {grid_fh})
        t_fh {str} -- table file handle (default: {table_fh})
        frmt {str} -- table format (default: {t_format})
        hdf5_pth {str} -- table hdf5 data path (default: {hdf5_pth})

    Returns:
        np.ndarray, astropy.table.Table -- grid, table
    '''
    return np.load(g_fh).T, Table.read(t_fh, format=frmt, path=hdf5_pth)

# =============================================================================
# =============================================================================


def regions_mpi(rank=0, size=1, r_0=r_start, r_1=r_stop, r_step=1,
                r_scale=annulus_scale, deg_0=-np.pi, deg_1=np.pi,
                deg_step=annulus_phi_step):
    '''

    [description]

    Keyword Arguments:
        rank {number} -- [description] (default: {0})
        size {number} -- [description] (default: {1})
        r_0 {number} -- [description] (default: {10})
        r_1 {number} -- [description] (default: {275})
        r_step {number} -- [description] (default: {1})
        r_scale {number} -- [description] (default: {0.1})
        deg_0 {number} -- [description] (default: {-180.0})
        deg_1 {number} -- [description] (default: {180.0})
        deg_step {number} -- [description] (default: {1.0})

    Returns:
        [type] -- [description]
    '''
    worklist = []
    t_rads = (r_1 - r_0) / r_step
    my_bit = int(t_rads / size)
    my_starting_point = int(my_bit * rank)
    for radius in range(my_starting_point, my_bit, r_step):

        region = radius * r_scale
        r_in = radius - region
        r_out = radius + region
        #r_in = radius
        #r_out = radius + r_step

        for annular_segment in np.linspace(
                deg_0,
                deg_1,
                deg_step,
                dtype=np.float16):

            worklist.append([r_in, r_out, deg_0, deg_1])

    return worklist


def record_table():
    # make an output table
    output_colums = ['radius', 'r0', 'r1', 'halo', 'xbox_cut',
                     'xbox_min', 'xbox_mean', 'xbox_max',
                     'run_length', 'sat_purity', 'n_stars',
                     'n_boxes', 'name0', 'name1', 'name2']
    output_dtyps = ['f', 'f', 'f', 'S10', 'f',
                    'f', 'f', 'f',
                    'i', 'f', 'i', 'i',
                    'f', 'f', 'f']
    return Table(names=output_colums, dtype=output_dtyps)


def fix_rslice(grid, r=4):
    phi_slice = np.zeros((601, 601, 1))
    grid = np.concatenate((grid, phi_slice), axis=2)

    center = 300
    grid = grid
    for i in range(grid.shape[0]):
        for q in range(grid.shape[1]):
            value = np.sqrt(
                np.square(i - center) + np.square(q - center))
            try:
                grid[i, q, r] = value
                grid[i, q, 5] = np.arctan2(i - 300, q - 300)
            except IndexError as e:
                pass

    grid[:, :, 5] = np.rad2deg(grid[:, :, 5], grid[:, :, 5]) + 180.0
    return grid


def mpi_run():

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    name = MPI.Get_processor_name()

    grid, table = load_halo

    for r_in, r_out, deg_0, deg_1 in regions_mpi(rank=0, size=1):

        grid_idx = np.nonzero(
            np.logical_and(
                grid[:, :, 4] >= r_in,
                grid[:, :, 4] < r_out))

        mu = grid[:, :, 0][grid_idx].mean()
        xbox_min = mu * 0.0005

        one_before = False
        run_length = 0

        min_xbox = []
        max_xbox = []
        mean_xbox = []

        alims = np.nonzero(
            np.logical_and(
                np.logical_and(
                    grid[:, :, 5] >= deg0,
                    grid[:, :, 5] < deg1),
                np.logical_and(
                    grid[:, :, 4] >= r_strt,
                    grid[:, :, 4] < r_stop)))

        xbox = ((grid[:, :, 0][alims] - mu) / mu)

        if xbox.mean() > xbox_min:

            if not one_before:

                run_length = 0

            run_length += 1
            one_before = True
            min_xbox.append(xbox.min())
            mean_xbox.append(xbox.mean())
            max_xbox.append(xbox.max())

        else:

            if one_before:

                row = [r,
                       r_strt,
                       r_stop,
                       grids[grid_number].split('_')[0],
                       xbox_min,
                       round(np.asarray(min_xbox).min(), 4),
                       round(np.asarray(mean_xbox).mean(), 4),
                       round(np.asarray(max_xbox).max(), 4),
                       run_length,
                       0.0,
                       grid[:, :, 0][alims].sum(),
                       len(alims[0]),
                       0.0,
                       0.0,
                       0.0]

                record_table.add_row(row)

            one_before = False
            run_length = 0
        m0 = '\n' + r + ' Kpc   '
        m1 = str(deg0) + ' deg  '
        m2 = '  '
        m3 = str(one_before) + ' '
        m4 = str(run_length) + ' '
        sys.stdout.write(m0 + m1 + m2 + m3 + m4)
        sys.stdout.flush()
