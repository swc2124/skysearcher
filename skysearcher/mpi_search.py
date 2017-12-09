from skysearch_lib import np
from skysearch_lib import os
from skysearch_lib import pause
from skysearch_lib import rank
from skysearch_lib import size
from skysearch_lib import name
from skysearch_lib import grid_list
from skysearch_lib import kpc_to_arcmin
from skysearch_lib import record_table
from skysearch_lib import radii
from skysearch_lib import save_record_table
from skysearch_lib import satid_setup
from skysearch_lib import load_grid
from skysearch_lib import time
from skysearch_lib import mu_idx
from skysearch_lib import get_annuli
from skysearch_lib import get_idx
from skysearch_lib import get_xbox
from skysearch_lib import new_sat_stars
from skysearch_lib import count_strs
from skysearch_lib import dom_satid


from skysearch_lib import XBOX_CUT
from skysearch_lib import N_SKIPS
from skysearch_lib import MIN_LOG_NSTARS
from skysearch_lib import MIN_N_SEGMENTS
from skysearch_lib import SAVE_INTERVAL

from skysearch_lib import stdout


#  Sat mass dict (m_book).
#  Sat acc_time (t_book).
m_book = np.load('../data/satmass_array.npy')
t_book = np.load('../data/satage_array.npy')

#  Misc values.
domsat_id = None
n_skips = 0
log_n_stars_max = [0]

#  Stall processes.
pause()
stdout.write(' '.join([
    'rank', 
    str(rank), 
    str(name), 
    str(size)]))
stdout.flush()


#  Load a list of halo names & grid file paths.
pause()
grids = grid_list()

# Read Mpc distance (d_mpc) from grid files.
dmpc_str = grids[0][1].split(os.path.sep)[-1].split('_')[1]
d_mpc = float(dmpc_str.replace('Mpc', ''))

# Make the mod for the distance.
mod = kpc_to_arcmin(d_mpc=d_mpc)

#  Load the record table (r_table).
pause()
r_table = record_table()

#  Load list of radii (radii).
#  i.e radii[x] = (r, r_start, r_stop)
pause()
radii = radii()

#  Designate work for rank (work_index).
pause()
work_index = range(rank, len(radii), size)
stdout.write('rank ' + str(rank) +  ' [ idx ] [' +
             ' '.join([str(i) for i in work_index]) + ']\n')
stdout.flush()

save_record_table(_table=r_table)

try:

    # Total run time(tic).
    tic = time()

    for halo, grid_fh in grids:

        # Get a list of satids and a table for counting
        # satids per region (satid_list) (satid_table).
        satid_list, satid_table = satid_setup(halo)

        # Get the data grid for this halo (grid).
        grid = load_grid(grid_fh)
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
            mu, r_idx = mu_idx(grid, r_start, r_stop)

            # The number of boxes in the annulus (log_n_boxes_in_ann).
            # log_n_boxes_in_ann = np.log10(len(np.nonzero(r_idx)[0]))

            # The number of stars in this annulus (log_n_stars_in_ann)
            # log_n_stars_in_ann = np.log10(grid[:, :, 0][r_idx].sum())

            # Boolean value to indicate the existence of
            # a immediately previous accepted segment(one_before).
            one_before = False

            # Load array of annuli segments and
            # step value (annuli) & (annuli_step).
            # i.e. annuli_step = float
            # i.e. annuli[x] = (deg_0, deg_1)
            annuli, annuli_step = get_annuli(r)

            # Number of empty segments (n_mt_seg).
            n_mt_seg = 0

            # Number of increasing segments (n_seg_increase).
            # Number of decreasing segments (n_seg_decrease).
            n_seg_increase = 0
            n_seg_decrease = 0
            last_nstars = 0

            # Step through the annulus (-pi, pi].
            # _deg0 = starting point.
            # _deg1 = ending point
            for _deg0, _deg1 in annuli:

                # The grid index for grid spaces within
                # this segment (idx).
                idx = get_idx(grid, _deg0, _deg1, r_idx)

                # The number of grid boxes this segment covers (n_boxes_tot).
                n_boxes_in_seg = len(idx[0])

                # If there are none, then continue.
                if not n_boxes_in_seg:
                    '''
                    stdout.write(
                        'rank ' 
                        + str(rank) 
                        +  ' [ EMPTY SEGMENT ] ' 
                        + halo 
                        + ' radius: ' + str(r) + ' Kpc'
                        + ' phi: ' + str(_deg0)
                        + '\n')
                    stdout.flush()
                    '''
                    n_mt_seg += 1
                    continue

                # The array of grid spaces from this segment's
                # contrast density value (xbox).
                xbox = get_xbox(grid, idx, mu)

                # Number of stars here (n_stars_here).
                n_stars_here = grid[:, :, 0][idx].sum()

                if n_stars_here >= last_nstars:
                    n_seg_increase += 1
                else:
                    n_seg_decrease += 1

                last_nstars = n_stars_here

                # ------------------------------------------------------
                #        Does this segment qualify to be a feature?
                # ------------------------------------------------------
                # Is the local min above xbox_cut?

                # If yes:
                if (
                    xbox.min() >= XBOX_CUT
                    and n_seg_increase >= n_seg_decrease
                    and np.log10(n_stars_here) >= MIN_LOG_NSTARS):
                    
                    # If this is a new feature...
                    if not one_before:

                        # open a list for the max values
                        # for this feature (xbmax).
                        xbmax = []

                        # open a list for max Log10(n_stars) (n_stars_max).
                        log_n_stars_max = []

                        # Start a dictionary for counting stars
                        # per satid (sat_stars).
                        sat_stars = new_sat_stars(satid_list)

                        # Set grid space count (n_stars).
                        n_boxes = 0

                        # Set star count (n_stars).
                        n_stars = 0

                        # Set the run length (n_segments) to 0.
                        n_segments = 0

                        # Remember the starting angular
                        # value (starting_deg).
                        starting_deg = _deg0

                        # Set allowed segment skips (n_skips).
                        n_skips = N_SKIPS

                    # Do this for every accepted segment:

                    # Put all region info into a list (r_info).
                    r_info = [r_start, r_stop, _deg0, _deg1]
                    
                    # Add this features min, mean and max
                    # to the existing lists. 
                    xbmax.append(xbox.max())

                    # Count stars per sat.
                    sat_stars = count_strs(
                        sat_stars,
                        r_info,
                        satid_table)

                    # Add boxes to total boxes.
                    n_boxes += n_boxes_in_seg

                    # Add stars to total stars for feature.
                    n_stars += n_stars_here

                    # Add n_stars to list for later.
                    log_n_stars_max.append(np.log10(n_stars_here))

                    # Increase run length (run_length) by 1.
                    n_segments += 1

                    # Remember that there is a feature
                    # currently being processed (one_before).
                    one_before = True


                # If no:
                else:

                    # Use allowed segment skips:
                    if n_skips >= 1:
                        n_skips -= 1
                        n_segments += 1
                        one_before = True

                    else:
                        
                        # If this is the end of a segment:
                        if (
                            one_before
                            and n_segments >= MIN_N_SEGMENTS):

                            # The feature's angular extent (angular_extent).
                            angular_extent = _deg0 - starting_deg

                            # The dominate satellite number (domsat_id).
                            # domsats's % of all stars (domsat_purity).
                            domsat_id, domsat_purity = dom_satid(sat_stars)

                            # An integer value for each halo (halo_num).
                            halo_num = int(halo[-2:])

                            # Mass of parent satellite (mass).
                            mass = m_book[domsat_id - 1]

                            # Accretion time of parent satellite (atime).
                            atime = t_book[domsat_id - 1]

                            # A new row for the r_table (row).
                            # Each feature is a row in the table.
                            row = [

                                # Halo.
                                halo_num,

                                # Annulus location values.
                                r, 
                                r_start, 
                                r_stop, 
                                annuli_step,
                                # n_mt_seg,

                                # Annulus content values.
                                # log_n_boxes_in_ann,
                                # log_n_stars_in_ann,
                                np.log10(mu),

                                # Feature content values.
                                max(xbmax),
                                max(log_n_stars_max),
                                domsat_purity,
                                domsat_id,
                                mass,
                                atime,

                                # Feature location values.
                                starting_deg,
                                # _deg0,
                                angular_extent,
                                # n_segments,
                                # n_boxes,

                                # rank
                                ]

                            r_table.add_row(row)

                            # Save point if added a new row.
                            if len(r_table) % SAVE_INTERVAL == 0:
                                save_record_table(_table=r_table)

                        # Clean up.
                        halo_num = None
                        mass = None
                        atime = None
                        angular_extent = 0
                        domsat_id = None
                        domsat_purity = None

                        # Reset grid space count (n_stars).
                        n_boxes = 0

                        # Reset n_stars to 0.
                        n_stars = 0

                        # Reset run length to 0.
                        n_segments = 0

                        # Remember that there's no previous segment
                        one_before = False

                        # Reset starting_deg.
                        starting_deg = None

                        n_seg_increase = 0
                        n_seg_decrease = 0
                        last_nstars = 0

            # Terminal progress message.
            line = ('rank ' + str(rank) +   ' : halo' + halo[-2:] + ' - ' +
                    str(r) + ' Kpc - ' +
                    str(round((time() - a_tic), 1)) + ' secs')
            stdout.write(line + '\n')
            stdout.flush()

        # Halo save point.
        save_record_table(_table=r_table)

    # Final save point.
    save_record_table(_table=r_table)

    # Exit message.
    msg = ('rank ' + str(rank) +    ' [ FINISHED ] [ ' +
           str(round((time() - tic) / 60.0, 1)) + ' minutes ]\n')
    stdout.write(msg)
    stdout.flush()
    exit(0)

except KeyboardInterrupt as e:
    print(e)
    exit(0)
