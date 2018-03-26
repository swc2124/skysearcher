:ref:`genindex`

skysearcher-newcfg
------------------
.. code-block:: bash

	$USER@machine: skysearcher-newcfg

	[ PATH ]
	data_dir = /lib/python3.6/site-packages/skysearcher/data
	plot_dir = /lib/python3.6/site-packages/skysearcher/data/plots
	table_dir = /lib/python3.6/site-packages/skysearcher/data/tables
	mpi_table_dir = /lib/python3.6/site-packages/skysearcher/data/tables/groupfinder/mpi
	grid_dir = /lib/python3.6/site-packages/skysearcher/data/grids
	grid_file_designator = grid
	grid_ext = npy
	table_file_designator = table
	table_format = hdf5
	table_hdf5_path = data
	table_ext = hdf5
	 
	[ Search Extent ]
	r_start = 5
	r_stop = 285
	r_step = 1
	annulus_scale = 0.05
	annulus_phi_step = 720
	 
	[ Accept Reject ]
	xbox_cut = 0.1
	min_log_nstars = 2.00
	min_n_segments = 4
	n_skips = 2
	 
	[ Run Time ]
	save_interval = 2
	 
	[ Data ]
	d_mpc = 4.0

:ref:`genindex`
