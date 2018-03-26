:ref:`genindex`

Making a new configuration file
-------------------------------
Once the configuration file is written, this is what you should have.  It will be in the skysearcher *root* directory named ``rc.cfg``.

.. code-block:: cfg

	[ PATH ]
	data_dir = ~/$USER/skysearcher/data
	plot_dir = ~/$USER/skysearcher/data/plots
	table_dir = ~/$USER/skysearcher/data/tables
	mpi_table_dir = ~/$USER/skysearcher/data/tables/groupfinder/mpi
	grid_dir = ~/$USER/skysearcher/data/grids
	
	grid_file_designator = grid
	table_file_designator = table
	table_format = hdf5
	table_hdf5_path = data
	table_ext = .hdf5
	grid_ext = npy
	 
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
