.. sectionauthor:: Sol Courtney <swc2124@columbia.edu>
.. codeauthor:: Sol Courtney <swc2124@columbia.edu>
.. py:module:: skysearcher

:ref:`genindex`

========================
**The search algorithm**
========================
This section will outline skysearcher parallel search algorithm step by step.

For each halo file
------------------

.. code-block:: python3
	:caption: The first loop introduces the halo grid file.  The grid file is a :py:mod:`numpy` array.
	:name: mpi_search.py

	for halo, grid_fh in grids:

		# Get a list of satids and a table for counting
		# satids per region (satid_list) (satid_table).
		satid_list, satid_table = satid_setup(halo)

		# Get the data grid for this halo (grid).
		grid = load_grid(grid_fh)

		# Adjust nstars for units and distance.
		grid[:, :, 1] /= mod

		# <commented out>
		# Kernel function
		# grid[:, :, 1] = gaussian_filter(grid[:,:,0], sigma=0.8, mode="constant", cval=0, order=0)

	""" End of first loop """

+ :term:`satid_list`, :term:`satid_table` = :func:`skysearcher.skysearch_lib.satid_setup`

+ :term:`grid` = :func:`skysearcher.skysearch_lib.load_grid`

+ :term:`grid` [:, :, 1] /= :term:`mod`

.. note:: commented out
	
	+ :term:`grid` [:, :, 1] = :func:`scipy.ndimage.filters.gaussian_filter`

For each radius assigned to this process
----------------------------------------
The second for loop concerns the radial selection.  Each process is assigned a list of radii for which they are responsible.  The second for loop iterates over this list.

.. code-block:: python3
	:lineno-start: 10

	""" End of first loop """
		
		# Step through the radii (r_start, r_stop).
		for job_id in work_index:

			# Start time (a_tic).
			a_tic = time()

			# r = radius in Kpc
			# r_start = starting radius (always < r)
			# r_stop = ending radius (always > r)
			r, r_start, r_stop = _radii[job_id]

			# This is so we dont need to index the whole table every
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

.. code-block:: python3

	# Step through the annulus (-pi, pi].
	# _deg0 = starting point.
	# _deg1 = ending point
	for _deg0, _deg1 in annuli:

		# Get the mu for this sub-annulus-section (mu).
		mu, r_idx2 = mu_idx2(grid, r_idx, _deg0, _deg1)

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
		xbox = get_xbox(grid, idx, mu)

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
		if (xbox.min() >= XBOX_CUT 
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
					if one_before and n_skips:
						n_skips -= 1
						n_segments += 1
						one_before = True

					else:

						# If this is the end of a segment:
						if (one_before and n_segments >= MIN_N_SEGMENTS):

							# The feature"s angular extent
							# (angular_extent).
							angular_extent = _deg0 - starting_deg

							# The dominate satellite number (domsat_id).
							# domsats"s % of all stars (domsat_purity).
							_va = dom_satid(sat_stars)
							domsat_id, domsat_purity, standout, nsats = _va

							# An integer value for each halo (halo_num).
							halo_num = int(halo[-2:])

							# Mass of parent satellite (mass).
							mass = m_book[domsat_id]

							# Accretion time of parent satellite (atime).
							atime = t_book[domsat_id]

							# Circle factor (jcirc).
							jcirc = j_book[domsat_id]

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
								standout,
								nsats,
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

:ref:`genindex`
