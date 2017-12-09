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
    import configparser
    import os

    # <setup>
    # Make a new configuration file parser object (config).
    config = configparser.ConfigParser()

    #
    # PATH
    # Get a string value for the path to the
    # current directory (curdir_path).
    curdir_path = os.path.abspath(os.path.curdir)

    # Get a string value for:
    # 1.) path to parent directory (pardir_path).
    # 2.) name of current directory (curdir_name).
    pardir_path, curdir_name = os.path.split(curdir_path)

    # Get a string value for parent directory name (pardir_name).
    pardir_name = os.path.split(pardir_path)[1]

    # Name the directory PATH for all data folders (data_dir).
    data_dir = os.path.join(pardir_path, 'data')
    table_dir = os.path.join(data_dir, 'tables')
    table_format = 'hdf5'
    config['PATH'] = {
        'data_dir': data_dir,
        'plot_dir': os.path.join(data_dir, 'plots'),
        'table_dir': table_dir,
        'mpi_table_dir': os.path.join(table_dir, 'groupfinder', 'mpi'),
        'grid_dir': os.path.join(data_dir, 'grids'),
        'grid_file_designator': 'grid',
        'grid_ext': 'npy',
        'table_file_designator': 'table',
        'table_format': table_format,
        'table_hdf5_path': 'data',
        'table_ext': os.path.extsep + table_format
        }

    config['Search Extent'] = {
        'r_start': '5',
        'r_stop': '285',
        'r_step': '1',
        'annulus_scale': '0.05',
        'annulus_phi_step': '720'
        }

    config['Accept Reject'] = {
        'xbox_cut': '0.1',
        'min_log_nstars': '2.00',
        'min_n_segments': '4',
        'n_skips': '2'
        }

    config['Run Time'] = {
        'save_interval': '2'
        }

    config['Data'] = {
        'd_mpc': '4.0'
        }
    
    if not rc_fh:
        # Make a new file name (file handle or fh).
        new_fh = 'rc' + os.path.extsep + 'cfg'

        # Join file name with the full path to current directory.
        rc_fh = os.path.join(curdir_path, new_fh)

    # Open file object (configfile) to write to.
    with open(rc_fh, 'w') as configfile:
        config.write(configfile)
    for sec in config.sections():
        print('[', sec, ']')
        for opt in config.options(sec):
            print(opt, '=', config.get(sec, opt))
        print(' ')
    return rc_fh

if __name__ == '__main__':
    new_rc()