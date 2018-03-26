# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


"""This is the main parallel skysearcher program.

Example
-------
    >>> mpiexec -nproc -machinefile python mpi_search.py

    -nproc {int} -- number of processes
    -machinefile {file} -- list of hosts each on a new line
"""

if __name__ == "__main__":
    from skysearch_lib import *
else:
    from .skysearch_lib import *

__all__ = ["main"]


def main():
    """This is the main parallel skysearcher program.

    Example
    -------
        >>> mpiexec -n <nproc> -machinefile <mf> python mpi_search.py

        -nproc {int} -- number of processes
        -machinefile {file} -- list of hosts each on a new line
    """

    # Load statice satillite data arrays.
    # Sat mass (m_book).
    # Sat acc_time (t_book).
    # Sat jcirc (j_book).
    m_book = np.load(os.path.join(DATA_DIR, "satmass_array.npy"))
    t_book = np.load(os.path.join(DATA_DIR, "satage_array.npy"))
    j_book = np.load(os.path.join(DATA_DIR, "satj_array.npy"))

    # Misc values.
    domsat_id = None
    n_skips = 0
    log_n_stars_max = [0]
    region_dict = {}

    # Stall processes.
    if MPI_RANK == 0:
        clear_tables(_target_dir=MPI_TABLE_DIR)


    # Print initial banner message.
    COMM.Barrier()
    pause()

    # Load a list of halo names & grid file paths.
    COMM.Barrier()
    pause()

    # Load file handles for numpy data arrays.
    grids = grid_list()

    # Read Mpc distance (d_mpc) from grid files.
    dmpc_str = grids[0][1].split(os.path.sep)[-1].split("_")[1]
    d_mpc = float(dmpc_str.replace("Mpc", ""))

    # Make the mod for the distance.
    mod = kpc_to_arcmin(d_mpc=d_mpc)

    # Load the record table (r_table).
    COMM.Barrier()
    pause()
    r_table = record_table(_names=TABLE_COLUMNS)

    # Load list of radii (radii).
    # i.e radii[x] = (r, r_start, r_stop)
    COMM.Barrier()
    pause()
    _radii = radii()

    # Designate work for MPI_RANK (work_index).
    COMM.Barrier()
    pause()
    work_index = range(MPI_RANK, len(_radii), MPI_SIZE)

    # Make print statements.
    COMM.Barrier()
    pause()

    # Flush Stdout buffer.
    COMM.Barrier()
    pause()

    # Total run time(tic).
    tic = time()

    # Human readable time.
    now = datetime.datetime.now()
    

    # Dictionary for printable information(run_dict).
    run_dict = {
        "MPI_RANK": MPI_RANK,
        "MPI_PROC_NAME": MPI_PROC_NAME,
        "START_TIME": now.strftime("%I:%M:%S %p - %A %B %d %Y"),
        "d_mpc": d_mpc,
        "mod": mod
    }
    
    # Save record hdf5 table.
    COMM.Barrier()
    pause()
    save_record_table(_table=r_table)
    report(_type="save", _info=run_dict)

    for halo, grid_fh in grids:

        h_tic = time()

        run_dict["halo"] = halo

        # Get a list of satids and a table for counting
        # satids per region (satid_list) (satid_table).
        satid_list, satid_table = satid_setup(halo)

        # Get the data grid for this halo (grid).
        grid = load_grid(grid_fh)

        # Print "starting halo" message.
        report(_type="starting halo", _info=run_dict)

        # Adjust nstars for units and distance.
        grid[:, :, 1] /= mod

        # <commented out>
        # Kernel function
        # grid[:, :, 1] = gaussian_filter(grid[:,:,0],
        # sigma=0.8, mode="constant", cval=0, order=0)

        # Step through the radii (r_start, r_stop).
        for job_id in work_index:

            # Start time (a_tic).
            a_tic = time()

            # r = radius in Kpc
            # r_start = starting radius (always < r)
            # r_stop = ending radius (always > r)
            r, r_start, r_stop = _radii[job_id]

            run_dict["r"] = r
            run_dict["r_start"] = r_start
            run_dict["r_stop"] = r_stop

            # This is so we don't need to index the whole table every
            # loop of the following for loop (local_satid_table).
            local_satid_table = satid_table[np.logical_and(
                satid_table["Rads"] >= r_start,
                satid_table["Rads"] < r_stop)]

            # Get the mu for this annulus (mu).
            mu, r_idx = mu_idx(grid, r_start, r_stop)

            # <commented out>
            # The number of boxes in the annulus (log_n_boxes_in_ann).
            # log_n_boxes_in_ann = np.log10(len(np.nonzero(r_idx)[0]))

            # <commented out>
            # The number of stars in this annulus (log_n_stars_in_ann)
            # log_n_stars_in_ann = np.log10(grid[:, :, 1][r_idx].sum())

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

                # Get the mu and indices's for this sub-annulus-section
                # (mu, r_idx2).
                mu_2, r_idx2 = mu_idx2(grid, r_idx, _deg0, _deg1)

                # The grid index for grid spaces within
                # this segment (idx).
                idx = get_idx(grid, _deg0, _deg1, r_idx2)

                # The number of grid boxes this segment covers
                # (n_boxes_tot).
                n_boxes_in_seg = len(idx[0])

                # If there are none, then continue.
                if not n_boxes_in_seg:
                    """
                    STDOUT.write(
                        "rank "
                        + str(MPI_RANK)
                        +  " [ EMPTY SEGMENT ] "
                        + halo
                        + " radius: " + str(r) + " Kpc"
                        + " phi: " + str(_deg0)
                        + "\n")
                    STDOUT.flush()
                    """
                    n_mt_seg += 1
                    continue

                # The array of grid spaces from this segment"s
                # contrast density value (xbox).
                xbox = get_xbox(grid, idx, mu_2)

                # Number of stars here (n_stars_here).
                n_stars_here = grid[:, :, 1][idx].sum()

                # TODO
                # What is this doing?
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
                    # -----------------------------------
                    # Put all region info into a list (r_info).
                    r_info = [r_start, r_stop, _deg0, _deg1]

                    # Add this features min, mean and max
                    # to the existing lists.
                    xbmax.append(xbox.max())

                    # Count stars per sat.
                    sat_stars = count_strs(
                        sat_stars,
                        r_info,
                        local_satid_table)

                    # Add boxes to total boxes.
                    n_boxes += n_boxes_in_seg

                    # Add stars to total stars for feature.
                    n_stars += int(n_stars_here)

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
                    if one_before and n_skips:
                        n_skips -= 1
                        n_segments += 1
                        one_before = True

                    else:

                        # If this is the end of a segment:
                        if (one_before 
                            and n_segments >= MIN_N_SEGMENTS
                            and sum(sat_stars.values())):


                            region_dict = {
                                "MPI_RANK": MPI_RANK
                            }

                            region_dict["r_start"] = r_start
                            region_dict["r_stop"] = r_stop
                            region_dict["annuli_step"] = annuli_step
                            region_dict["log(mu_2)"] = np.log10(mu_2)
                            region_dict["max xbox"] = max(xbmax)
                            region_dict["max log(nstr)"] = max(log_n_stars_max)
                            region_dict["n_stars"] = int(n_stars)
                            region_dict["n_segments"] = n_segments
                            region_dict["n_seg_increase"] = n_seg_increase
                            region_dict["n_seg_decrease"] = n_seg_decrease

                            # The dominate satellite number (domsat_id).
                            # domsats"s % of all stars (domsat_purity).
                            region_dict = dom_satid(sat_stars, region_dict)


                            # The feature"s angular extent
                            # (angular_extent).
                            angular_extent = _deg0 - starting_deg
                            region_dict["starting deg"] = starting_deg
                            region_dict["angular extent"] = angular_extent

                            # An integer value for each halo (halo_num).
                            halo_num = int(halo[-2:])
                            # region_dict["halo_num"] = halo_num
                            region_dict["halo"] = halo

                            # Mass of parent satellite (mass).
                            mass = m_book[region_dict["domsat_id"]]
                            region_dict["mass"] = mass

                            # Accretion time of parent satellite (atime).
                            atime = t_book[region_dict["domsat_id"]]
                            region_dict["atime"] = atime

                            # Circle factor (jcirc).
                            jcirc = j_book[region_dict["domsat_id"]]
                            region_dict["jcirc"] = jcirc

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
                                np.log10(mu_2),

                                # Feature content values.
                                max(xbmax),
                                max(log_n_stars_max),
                                region_dict["domsat_purity"],
                                region_dict["domsat_id"],
                                region_dict["standout"],
                                region_dict["n_sats_end"],
                                mass,
                                atime,
                                jcirc,

                                # Feature location values.
                                starting_deg,
                                # _deg0,
                                angular_extent,
                                # n_segments,
                                # n_boxes,

                                # MPI values.
                                # MPI_RANK
                            ]

                            # Add the row to the table.
                            r_table.add_row(row)

                            # Save point if added a new row.
                            if len(r_table) % SAVE_INTERVAL == 0:
                                save_record_table(_table=r_table)
                                report(_type="save", _info=run_dict)


                            if MPI_RANK%3==0:
                                report(_type="new region", _info=region_dict)
                            del region_dict

                        # Clean up.
                        halo_num = None
                        mass = None
                        atime = None
                        angular_extent = 0

                        # Reset grid space count (n_stars).
                        n_boxes = 0

                        # Reset n_stars to 0.
                        n_stars = 0

                        # Reset run length to 0.
                        n_segments = 0

                        # Remember that there"s no previous segment
                        one_before = False

                        # Reset starting_deg.
                        starting_deg = None
                        n_seg_increase = 0
                        n_seg_decrease = 0
                        last_nstars = 0



            # Terminal progress message.
            run_dict["annulus time"] = round((time() - a_tic))
            report(_type='end annulus', _info=run_dict)

        # Halo save point.
        save_record_table(_table=r_table)
        report(_type="save", _info=run_dict)
        run_dict["halo time"] = round((time() - h_tic) / 60, 1)
        report(_type='end halo', _info=run_dict)


    # Final save point.
    save_record_table(_table=r_table)
    report(_type="save", _info=run_dict)

    # Exit message.
    run_dict["program time"] = round((time() - tic) / 60, 1)
    report(_type='exit', _info=run_dict)
    exit(0)


if __name__ == "__main__":
    main()
