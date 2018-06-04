
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from time import time

from mpi4py import MPI

import skysearch_lib as ss_lib

from skysearch_lib import np
from skysearch_lib import os
from skysearch_lib import sleep

stdout = os.sys.stdout

# C:\Users\swc21\GitHub\skysearcher\skysearcher\mpi_search.py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

#  Stall processes.
sleep(rank * 0.005)
stdout.write(' '.join(['rank', str(rank), str(name), str(size)]))
stdout.flush()

#  Sat mass dict (m_book).
m_book = np.load('../data/satmass_array.npy')
t_book = np.load('../data/satage_array.npy')
#

#  Start time (tic).
tic = time()

#  Misc values.
dom_sat = None
sat_fraction = 0.001

#  Load a list of halo names & grid file paths.
sleep(0.005)
grids = ss_lib.grid_list()

# Read Mpc distance (d_mpc) from grid files.
d_mpc = float(grids[0][1].split(os.path.sep)[-1].split('_')[1].replace('Mpc', ''))

# Make the mod for the distance.
mod =ss_lib.kpc_to_arcmin(d_mpc=d_mpc)
#  Load the record table (r_table).
sleep(0.005)
r_table = ss_lib.record_table()

#  File naming regime.
r_table.meta['fh'] = os.path.join(ss_lib.DATA_DIR,
                                  'tables',
                                  'groupfinder',
                                  'mpi',
                                  str(rank) + '.hdf5')
r_table.meta['xbox_min_fixed'] = ss_lib.XBOX_MIN_FIXED
r_table.meta['xbox_min_mu'] = ss_lib.XBOX_MIN_MU
#  Load list of radii (radii).
#  i.e radii[x] = (r, r_start, r_stop)
sleep(rank * 0.015)
radii = ss_lib.radii()

#  Designate work for rank (work_index).
sleep(0.005)
work_index = range(rank, len(radii), size)
stdout.write('rank ' + str(rank) +  ' [ idx ] [' +
             ' '.join([str(i) for i in work_index]) + ']\n')
stdout.flush()
n_jobs = len(work_index)

sleep(0.005)

def save_record_table(_table=r_table, _rnk=rank):
    '''Helper function to safely write table to disc.

    Arguments:
        _table {dict} -- This is the record_table
    '''
    try:
        # Write the table to disc.
        _table.write(_table.meta['fh'],
                     format='hdf5',
                     path='data',
                     overwrite=True,
                     append=True)
        stdout.write('rank ' + str(_rnk) + ' [ SAVED ]\n')
        stdout.flush()
    except IOError as e:
        # Try again in a second.
        stdout.write('[ ' + str(_rnk) + ' ] [ ERROR ] ' + e.args[0] + '\n')
        sleep(1)
        save_record_table(_table=_table)

try:
    counter = 0
    save = 20
    for halo, grid_fh in grids:

        #halo, grid_fh = grids[-2]

        stdout.write('rank ' + str(rank) +  ' -> ' + grid_fh.split('\\')[-1] + '\n')
        stdout.flush()

        # Get a list of satids and a table for counting
        # satids per region (satid_list) (satid_table).
        satid_list, satid_table = ss_lib.satid_setup(halo)

        # Plot a full polar projection of the halo
        # for reference.
        # ss_lib.plot_full_halo(halo)

        # make an empty plot to add segments to (data_plot).
        # data_plot = ss_lib.get_data_plot(halo)

        # Get the data grid for this halo (grid).
        grid = ss_lib.load_grid(grid_fh)
        stdout.write('rank ' + str(rank) +  ' [ LOADED ] ' + halo + '\n')
        stdout.flush()

        # Adjust nstars for units and distance.
        grid[:, :, 0] /= mod

        # Step through the radii (r_start, r_stop).
        for job_id, _i in enumerate(work_index):

            # Start time (a_tic).
            a_tic = time()

            # r = radius in Kpc
            # r_start = starting radius (always < r)
            # r_stop = ending radius (always > r)
            r, r_start, r_stop = radii[_i]

            # get the mu for this annulus (mu).
            mu, r_idx = ss_lib.mu_idx(grid, r_start, r_stop)

            # Boolean value to indicate the existence of
            # a immediately previous accepted segment(one_before).
            one_before = False
            skips = 4

            # An integer value for the number of consecutively
            # accepted regions (run_length).
            run_length = 0

            # Load array of annuli segments and
            # step value (annuli) & (annuli_step).
            # i.e. annuli_step = float
            # i.e. annuli[x] = (deg_0, deg_1)
            annuli, annuli_step = ss_lib.annuli(r)
            #EXP_n_steps = len(annuli)

            # Step through the annulus (-pi, pi].
            # _deg0 = starting point.
            # _deg1 = ending point
            for _deg0, _deg1 in annuli:

                #swath = np.pi / 4.0
                # get the mu for this annulus (mu).
                # mu, r_idx = ss_lib.mu_idx(grid, r_start, r_stop)
                #mu, r_idx = ss_lib.EXP_mu_idx(grid, _deg0 - swath, _deg1 + swath, r_start, r_stop)

                # The grid index for grid spaces within
                # this segment (idx).
                idx = ss_lib.get_idx(grid, _deg0, _deg1, r_idx)

                # EXPERIMENTAL XBOX INDEXING TEST
                #os.system('cls')
                #idx = ss_lib.EXPERIMENTAL_IDX(grid, _deg0, _deg1, r_idx, annuli_step, EXP_n_steps)

                if not len(idx[0]):
                    skips -= 1
                    if skips == 0:
                        one_before = False
                    continue

                # The array of grid spaces from this segment's
                # contrast density value (xbox).
                xbox = ss_lib.xbox(grid, idx, mu)
                # xbox = (grid[:, :, 0][idx] - mu) / mu

                # ------------------------------------------------------
                #        Does this segment qualify to be a feature?
                # ------------------------------------------------------
                # If it does,
                # do this:
                # Is the local mean above xbox_min?
                _xbpass = xbox.min() # + (2.0 * xbox.std())
                #print('ss_lib.XBOX_MIN :', ss_lib.XBOX_MIN)
                #print('_xmean          :', round(_xmean, 4))
                #print('RADIUS          :', r, 'Kpc')
                #sleep(.1)

                _xbox_mu = ss_lib.XBOX_MIN_MU / np.sqrt(mu)
                _xbox_fixed = ss_lib.XBOX_MIN_FIXED

                xbox_min = min([ss_lib.XBOX_MIN, _xbox_mu, _xbox_fixed])
                if _xbpass >= ss_lib.XBOX_MIN:

                    # Put all region info into a list (r_info).
                    r_info = [r_start, r_stop, _deg0, _deg1]

                    # Make a new feature.        <--
                    if not one_before or run_length > ss_lib.MAX_SEG_SIZE:

                        # open a list for the min, mean, max values
                        # for this feature (xbmin, xbmean, xbmax).
                        xbmin = []
                        xbmean = []
                        xbmax = []

                        # open a list for max Log10(n_stars)
                        n_stars_max = []

                        # Start a dictionary for counting stars
                        # per satid (sat_stars).
                        sat_stars = ss_lib.new_sat_stars(satid_list)

                        # Reset grid space count (n_stars).
                        n_boxes = 0

                        # Reset star count (n_stars).
                        n_stars = 0

                        # Reset the run length (run_length) to 0.
                        run_length = 0

                        # Remember the starting angular
                        # value (starting_deg).
                        starting_deg = _deg0

                    # Add this features min, mean and max
                    # to the existing lists.
                    xbmin.append(xbox.min())
                    xbmean.append(xbox.mean())
                    xbmax.append(xbox.max())

                    # Count stars per sat.
                    sat_stars = ss_lib.count_strs(
                        sat_stars,
                        r_info,
                        satid_table)

                    # Add boxes to total boxes.
                    n_boxes += len(idx[0])

                    # Add stars to total stars for feature.
                    _n_stars_here = grid[:, :, 0][idx].sum()
                    n_stars += _n_stars_here

                    # Add n_stars to list for later.
                    n_stars_max.append(_n_stars_here)

                    # Increase run length (run_length) by 1.
                    run_length += 1

                    # Remember that there is a feature
                    # currently being processed (one_before).
                    one_before = True

                # If it does not,
                # do this:
                else:
                    
                    if one_before:
                        # The feature's angular extent (angular_extent).
                        angular_extent = _deg0 - starting_deg
                        # If a feature is ending,
                        # do this:

                        # The dominate satellite number (dom_sat).
                        # dom_sats's % of all stars (n_stars).
                        dom_sat, sat_frc = ss_lib.dom_satid(sat_stars)

                        if (
                            dom_sat
                            and np.log10(n_stars) >= ss_lib.MIN_LOG_NSTARS
                            #and angular_extent <= 0.25 * np.pi
                            and run_length >= ss_lib.MIN_SEG_SIZE
                            #and run_length <= ss_lib.MAX_SEG_SIZE):
                            ):

                            # An integer value for each halo (halo_num).
                            halo_num = int(halo[-2:])

                            # Mass of parent satellite (mass).
                            mass = m_book[dom_sat - 1]

                            # Accretion time of parent satellite (atime).
                            atime = t_book[dom_sat - 1]

                            # A new row for the r_table (row).
                            # Each feature is a row in the table.
                            row = [
                                # Search region values.
                                r, r_start, r_stop, annuli_step,

                                # Halo
                                halo_num,

                                # Annulus values.
                                xbox_min, _xbox_fixed, _xbox_mu, np.log10(mu),
                                
                                # Feature values.
                                min(xbmin),
                                np.asarray(xbmean).mean(),
                                max(xbmax),
                                np.log10(max(n_stars_max)),
                                np.log10(n_stars), 

                                sat_frc, dom_sat, mass, atime,
                                
                                starting_deg, _deg0, angular_extent,
                                run_length, n_boxes,

                                # Process values.
                                rank
                            ]
                            r_table.add_row(row)

                            counter += 1
                            if save == counter:
                                save_record_table(_table=r_table)
                                counter = 0


                            # Write a log file for the stars in dom_sat.
                            # ss_lib.log_satstars(sat_stars, halo, r, dom_sat)

                            # Add the segment area to the data_plot & save.
                            # ---------------------------------------------
                            # Left most point (theta0).
                            # or: _deg1 - (0.5 * angular_extent)
                            # theta0 = starting_deg

                            # The angular extent (angular_extent).
                            # angular_extent = angular_extent

                            # The starting radius value (r_start).
                            # r_start = r_start

                            # The radial extent (r_extent).
                            # r_extent = r_stop - r_start

                            # A list of location information needed
                            # to plot (_pts).
                            # _pts = [theta0, r_extent, angular_extent, r_start]
                            # ss_lib.update_plot(data_plot, _pts, halo, dom_sat)

                            # Clean up.
                            halo_num = None
                            mass = None
                            del row

                        # Clean up.
                        angular_extent = 0
                        dom_sat = None
                        sat_frc = None

                        # Reset grid space count (n_stars).
                        n_boxes = 0

                        # Reset n_stars to 0.
                        n_stars = 0

                        # Reset run length to 0.
                        run_length = 0

                        # Remember that there's no previous segment
                        one_before = False

                        # Reset starting_deg.
                        starting_deg = None

            # Terminal progress message.
            line = ('rank ' + str(rank) +   ' : halo' + halo[-2:] + ' - ' +
                    str(r) + ' Kpc - ' +
                    str(round((time() - a_tic), 1)) + ' secs')
            stdout.write(line + '\n')
            stdout.flush()

            # Save table.
            # save_record_table(_table=r_table)
            # save = False
            # save_record_table(_table=r_table)
        #break

        # Save the table.
        save_record_table(_table=r_table)

    # Save table.
    save_record_table(_table=r_table)
    msg = ('rank ' + str(rank) +    ' [ FINISHED ] [ ' +
           str(round((time() - tic) / 60.0, 1)) + ' minutes ]\n')
    stdout.write(msg)
    stdout.flush()

except KeyboardInterrupt as e:
    print(e)
    exit(0)
