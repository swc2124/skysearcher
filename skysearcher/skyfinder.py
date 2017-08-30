'''[summary]

[description]

Variables:
    record_table {[type]} -- [description]
    for fh in grids: {[type]} -- [description]
'''
from __future__ import division, absolute_import, print_function

from sys import stdout

import numpy as np
import os

import skysearch_lib as ss_lib

dom_sat = None
sat_fraction = 1
# Load a list of halo names & grid file paths.
grids = ss_lib.grid_list()

# Load the record table (record_table).
record_table = ss_lib.record_table()
table_fh = os.path.join(ss_lib.table_dir, 'groupfinder', 'main2.hdf5')
record_table.write(table_fh, format='hdf5', path='data', overwrite=True)

# Load list of radii (radii).
# i.e radii[x] = (r, r_start, r_stop)
radii = ss_lib.radii()

# Load array of annuli segments and
# step value (annuli) & (annuli_step).
# i.e. annuli_step = float
# i.e. annuli[x] = (deg_0, deg_1)
annuli, annuli_step = ss_lib.annuli()

for halo, grid_fh in grids[3:]:

    # Get a list of satids and a table for counting
    # satids per region (satid_list) (satid_table).
    satid_list, satid_table = ss_lib.satid_setup(halo)

    # Plot a full polar projection of the halo
    # for reference.
    # ss_lib.plot_full_halo(halo)

    # make an empty plot to add segments to (data_plot).
    #data_plot = ss_lib.get_data_plot(halo)

    # Get the data grid for this halo (grid).
    grid = ss_lib.load_grid(grid_fh)

    # Step through the radii (r_start, r_stop).
    # i.e. r = radius
    # i.e. r_start = starting radius (always < r)
    # i.e. r_stop = ending radius (always > r)
    for r, r_start, r_stop in radii:

        # get the mu for this annulus (mu).
        mu, r_idx = ss_lib.mu_idx(grid, r_start, r_stop)

        # Set the minimum value a segment must have to
        # be considered for a feature (xbox_min).
        #xbox_min = ss_lib.xbox_min
        xbox_min = 5.0 / np.sqrt(mu)

        # Boolean value to indicate the existence of
        # a immediately previous accepted segment(one_before).
        one_before = False

        # An integer value for the number of consecutively
        # accepted regions (run_length).
        run_length = 0

        # Step through the annulus.
        # i.e. _deg0 = starting point.
        # i.e. _deg1 = ending point
        for _deg0, _deg1 in annuli:

            # The grid index for grid spaces within
            # this segment (seg_idx).
            idx = ss_lib.get_idx(grid, _deg0, _deg1, r_idx)

            # the array of grid spaces from this segment's
            # contrast density value (xbox).
            xbox = ss_lib.xbox(grid, idx, mu)
            #xbox = (grid[:, :, 0][idx] - mu) / mu

            # ------------------------------------------
            # Does this segment qualify to be a feature?
            # ------------------------------------------
            # If it does,
            # do this:
            if xbox.mean() > xbox_min:

                # Put all region info into a list (r_info).
                r_info = [r_start, r_stop, _deg0, _deg1]

                # Make a new feature.        <--
                if not one_before:

                    # open a list for the min, mean, max values
                    # for this feature (xbmin, xbmean, xbmax).
                    xbmin = []
                    xbmean = []
                    xbmax = []

                    # Start a dictionary for counting stars
                    # per satid (sat_stars).
                    sat_stars = ss_lib.new_sat_stars(satid_list)

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

                # Add stars to total stars for feature.
                n_stars += grid[:, :, 0][idx].sum()

                # Increase run length (run_length) by 1.
                run_length += 1

                # Remember that there is a feature
                # currently being processed (one_before).
                one_before = True

            # If it does not,
            # do this:
            else:

                # If a feature is ending,
                # do this:
                if one_before and run_length > ss_lib.min_seg_size:

                    # The feature's angular extent (angular_extent).
                    angular_extent = _deg0 - starting_deg

                    dom_sat, sat_fraction = ss_lib.dom_satid(sat_stars)

                    # A new row for the record_table (row).
                    # Each feature is a row in the table.
                    row = [r, r_start, r_stop,
                           halo, xbox_min,
                           np.asarray(xbmin).mean(),
                           np.asarray(xbmean).mean(),
                           np.asarray(xbmax).mean(),
                           angular_extent,
                           sat_fraction,
                           n_stars,
                           len(idx[0]),
                           dom_sat,
                           0.0, 0.0]
                    record_table.add_row(row)

                    # Add the segment area to the data_plot & save.
                    # ---------------------------------------------
                    # Left most point (theta0).
                    # or: _deg1 - (0.5 * angular_extent)
                    #theta0 = starting_deg

                    # The angular extent (angular_extent).
                    # angular_extent = angular_extent

                    # The starting radius value (r_start).
                    # r_start = r_start

                    # The radial extent (r_extent).
                    #r_extent = r_stop - r_start

                    # A list of location information needed
                    # to plot (_pts).
                    #_pts = [theta0, r_extent, angular_extent, r_start]
                    #ss_lib.update_plot(data_plot, _pts, halo, dom_sat)

                # Remember that there's no previous segment
                one_before = False

                # Reset run length to 0.
                run_length = 0

        stdout.write('\r' + halo + '     r : ' + str(r) +
                     ' Kpc      mu : ' + str(mu) + '    xboxmin : ' +
                     str(xbox_min) + '     domsat : ' + str(dom_sat) +
                     '    %  : ' + str(round(sat_fraction * 1e2, 2)))
        stdout.flush()

        #ss_lib.save_plot(data_plot, halo)
        if r in range(0, 300, 10):
            record_table.write(
                table_fh,
                format='hdf5',
                path='data',
                overwrite=True)
        #print('saved', halo, 'table')
