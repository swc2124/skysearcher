'''
===============================================================================
skysearcher v0.1.0
August 2017
Sol W. Courtney
swc2124@columbia.edu
Columbia University NYC, NY
===============================================================================
python skysearcher/new_cfg.py

This script contains a single function <new_rc()>.  It is to be called
in the event that no configuration file exists in the main skysearcher
directory <skysearcher/skysearcher>.  After running this function, there
will be a new configuration file titled "rc.cfg."  Skysearcher will read
this file.  The file can be manually edited as text prior to running
skysearcher.

to read this file:
config = ConfigParser.RawConfigParser()
config.read('rc.cfg')

ConfigParser Documentation:
https://docs.python.org/2/library/configparser.html
'''


def new_rc(rc_fh=None):
    '''
    create and write a new rc configuration file to
    the provided file handle

    Keyword Arguments:
        rc_fh {str} -- full path and name of rc file (default: {None})

    Returns:
        rc file handle str.
    '''
    # <imports>
    import ConfigParser
    import os

    # <setup>
    # Make a new configuration file parser object (config).
    config = ConfigParser.ConfigParser()

    #
    # <select data parameters ("Data")>
    # =========================================================================
    # dmpc = '4.0'
    # filter_type = 'h158'
    # halo ='halo02'
    # file_type = 'grid'
    # ext = 'npy'
    config.add_section('Data')

    # Set the data cut level (data_cut). Example data[::cut]
    config.set('Data', 'data_cut', '1')

    # Name the halo to be loaded (halo).
    config.set('Data', 'halo', 'halo02')

    # Set the distance in mega parsecs (d_mpc).
    config.set('Data', 'd_mpc', '4.0')

    # Set the WFIRST filter type (filter_type).
    config.set('Data', 'filter_type', 'h158')

    #
    # <select or set PATH values ("PATH")>
    # =========================================================================
    config.add_section('PATH')

    # Get a string value for the path to the
    # current directory (curdir_path).
    curdir_path = os.path.abspath(os.path.curdir)
    config.set('PATH', 'curdir_path', curdir_path)

    # Get a string value for:
    # 1.) path to parent directory (pardir_path).
    # 2.) name of current directory (curdir_name).
    pardir_path, curdir_name = os.path.split(curdir_path)
    config.set('PATH', 'curdir_name', curdir_name)
    config.set('PATH', 'pardir_path', pardir_path)

    # Get a string value for parent directory name (pardir_name).
    pardir_name = os.path.split(pardir_path)[1]
    config.set('PATH', 'pardir_name', pardir_name)

    # Name the directory PATH for all data folders (data_dir).
    data_dir = os.path.join(pardir_path, 'data')
    config.set('PATH', 'data_dir', data_dir)

    # Name the directory PATH for all plots (plot_dir).
    plot_dir = os.path.join(data_dir, 'plots')
    config.set('PATH', 'plot_dir', plot_dir)

    # Name the directory PATH for all tables (table_dir).
    table_dir = os.path.join(data_dir, 'tables')
    config.set('PATH', 'table_dir', table_dir)

    # Name the directory PATH for all grids (grid_dir).
    grid_dir = os.path.join(data_dir, 'grids')
    config.set('PATH', 'grid_dir', grid_dir)

    # Set the base file handle prefix (fh_prefix).
    # Example: 'halo02_4.0_h158' i.e halo_dmpc_filter_type
    fh_prefix = '_'.join([
        config.get('Data', 'halo'),
        config.get('Data', 'd_mpc'),
        config.get('Data', 'filter_type')])
    config.set('PATH', 'fh_prefix', fh_prefix)

    # <<GRID FILES>>
    # Grid file designator (grid_file_designator).
    config.set('PATH', 'grid_file_designator', 'grid')

    # Grid file extension (grid_ext).
    config.set('PATH', 'grid_ext', os.path.extsep + 'npy')

    # <<TABLE FILES>>
    # Table file designator (table_file_designator).
    config.set('PATH', 'table_file_designator', 'table')

    # Table file type (table_format).
    table_format ='hdf5'
    config.set('PATH', 'table_format', table_format)

    # Table hdf5 path (table_hdf5_path).
    config.set('PATH', 'table_hdf5_path', 'data')

    # Table file extension (table_ext).
    config.set('PATH', 'table_ext', os.path.extsep + table_format)



    #
    # <select distance units ("Program Units")>
    # =========================================================================
    config.add_section('Program Units')

    # Set the units to be used for the halo's radius (units_radius).
    config.set('Program Units', 'units_radius', 'kpc')

    # Set the units to be used for rotating about the
    # annulus (units_annulus).
    config.set('Program Units', 'units_annulus', 'degree')

    #
    # <select search extent ("Search Extent")>
    # =========================================================================
    config.add_section('Search Extent')

    # Set the beginning radius (r_start).
    config.set('Search Extent', 'r_start', 10)

    # Set the ending radius (r_stop).
    config.set('Search Extent', 'r_stop', 275)

    # Set the percentage of radius to be used as the extent of the
    # annulus (annulus_scale).
    # a0, a1 = radius + and - (radius * annulus_scale) / 2.0
    config.set('Search Extent', 'annulus_scale', 0.1)

    # Set the step length for each step around the
    # annulus (annulus_phi_step).
    config.set('Search Extent', 'annulus_phi_step', 1)

    # Set the minimum value for xbox (contrast density) to be considered
    # for the search (xbox_min).
    config.set('Search Extent', 'xbox_min', 10.0)

    #
    # <select accept-reject criterion ("Accept Reject")>
    # =========================================================================
    config.add_section('Accept Reject')

    # Set the minimum number of grid spaces (boxes) for
    # features (min_boxes).
    # Prospect features must have at least this many boxes to be
    # accepted as a feature.
    config.set('Accept Reject', 'min_boxes', 25)

    #
    # <select plots ("Plots")>
    # =========================================================================
    config.add_section('Plots')

    # Plot a full projected image of the halo using each stars px, py
    # (full_projection)? True/False
    config.set('Plots', 'full_projection', True)

    # Plot a heat-map for the halo grid (heatmap)?
    # True/False
    config.set('Plots', 'heatmap', True)

    # Plot all regions together on a single plot (True), or plot them
    # individually (False)(regions_together)?
    # True/False
    config.set('Plots', 'regions_together', True)

    #
    # <write rc file to disc>
    # =========================================================================
    if not rc_fh:

        # Make a new file name (file handle or fh).
        new_fh = 'rc' + os.path.extsep + 'cfg'

        # Join file name with the full path to current directory.
        rc_fh = os.path.join(curdir_path, new_fh)

    # Open file object (configfile) to write to.
    with open(rc_fh, 'wb') as configfile:
        config.write(configfile)

    return rc_fh
