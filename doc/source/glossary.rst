.. sectionauthor:: Sol Courtney <swc2124@columbia.edu>
.. codeauthor:: Sol Courtney <swc2124@columbia.edu>

:ref:`genindex`

============
**Glossary**
============

.. <BASIC TERMS> 2018-03-24 ##################################################
.. glossary::

	Arc-min : Units

		One minute of arc.
	
	Kpc : Units

		:math:`10^3` parsecs or Kilo-parsec.

	Mpc : Units

		:math:`10^6` parsecs or Mega-parsec.

	Radians : Units

		A unit of angle, equal to an angle at the center of a circle whose arc is equal in length to the radius.

.. <CONFIGURATION> 2018-03-22 ################################################
.. glossary::

	.. <CONFIGURATION> 2018-03-22 21:20:24
	DATA_DIR : Configfile-PATH
		
		Path to skysearcher's input and output data directory.  This should be an absolute path.

		Example: ~/$USER/lib/python3.6/site-packages/skysearcher/data/plots
		
		Set in :ref:`config-PATH`.

	PLOT_DIR : Configfile-PATH
		
		Path to plot output directory.  This should be an absolute path.

		Example: ~/$USER/lib/python3.6/site-packages/skysearcher/data/plots
			
		Set in :ref:`config-PATH`.

	TABLE_DIR : Configfile-PATH
		
		Path to hdf5 table input and output directory.  This should be an absolute path.

		Example: ~/$USER/lib/python3.6/site-packages/skysearcher/data/tables
			
		Set in :ref:`config-PATH`.

	MPI_TABLE_DIR : Configfile-PATH
		
		Path to mpi hdf5 table directory within :term:`TABLE_DIR`.  This should be an absolute path.

		Example: ~/$USER/lib/python3.6/site-packages/skysearcher/data/tables/groupfinder/mpi
			
		Set in :ref:`config-PATH`.

	GRID_DIR : Configfile-PATH
		
		Path to numpy array input files directory.  This should be an absolute path.

		Example: ~/$USER/lib/python3.6/site-packages/skysearcher/data/grids
			
		Set in :ref:`config-PATH`.

	TABLE_FORMAT : Configfile-PATH

		Table format for :py:mod:`astropy.tables.Table` (default: {"hdf5"})
			
		Set in :ref:`config-PATH`.

	TABLE_HDF5_PATH : Configfile-PATH
		
		Table path for hdf5 table (default: {"hdf5"})
			
		Set in :ref:`config-PATH`.

	TABLE_EXT : Configfile-PATH
		
		Table file-type extension (default: {"hdf5"})
			
		Set in :ref:`config-PATH`.

	GRID_EXT : Configfile-PATH
		
		Extension for array files. (default: {"npy"})
			
		Set in :ref:`config-PATH`.

	R_START : Configfile-Search-Extent
		
		Starting radius in units of :term:`Kpc` to start search from.
			
		Set in :ref:`config-Search Extent`.

	R_STOP : Configfile-Search-Extent
		
		Ending search radius in units of :term:`Kpc`.
			
		Set in :ref:`config-Search Extent`.

	R_STEP : Configfile-Search-Extent
		
		Step interval in units of :term:`Kpc`.
			
		Set in :ref:`config-Search Extent`.

	R_SCALE : Configfile-Search-Extent
		
		Percent to multiply radius by to get + and - :term:`Kpc` for :term:`annulus`.
			
		Set in :ref:`config-Search Extent`.

	ANNULUS_PHI_STEP : Configfile-Search-Extent
		
		How many sections to divide each annulus by.
			
		Set in :ref:`config-Search Extent`.

	XBOX_CUT : Configfile-Accept-Reject
		
		Minimum contrast density a feature must have.
			
		Set in :ref:`config-Accept Reject`.

	MIN_LOG_NSTARS : Configfile-Accept-Reject
		
		Minimum log10 value for number of stars a feature must have.
			
		Set in :ref:`config-Accept Reject`.

	MIN_N_SEGMENTS : Configfile-Accept-Reject
		
		Minimum number of segments a feature must be
			
		Set in :ref:`config-Accept Reject`.

	N_SKIPS : Configfile-Accept-Reject
		
		Allowed number of skipped segments before ending feature.
			
		Set in :ref:`config-Accept Reject`.

	SAVE_INTERVAL : Configfile-Run-Time
		
		Number of features to hold between saves.
			
		Set in :ref:`config-Run Time`.

	D_MPC : Configfile-Data
		
		Distance to target halo in units of :term:`Mpc`.
			
		Set in :ref:`config-DATA`.
	.. <END - CONFIGURATION>

.. <SKYSEARCH_LIB> 2018-03-24 ################################################
.. glossary::

	.. <SKYSEARCH_LIB>
	COMM : skysearch_lib.py

		MPI COMM_WORLD communicator.
		
		:py:mod:`MPI.COMM_WORLD`

	MPI_RANK : skysearch_lib.py

		A unique integer value for each MPI process. 
		
		:py:func:`COMM.Get_rank()`

	MPI_SIZE : skysearch_lib.py

		Number of MPI processes in use.
		
		:py:func:`COMM.Get_size()`

	MPI_PROC_NAME : skysearch_lib.py

		Name of MPI process's host machine.
		
		:py:func:`MPI.Get_processor_name()`

	START_TIME : skysearch_lib.py
		
		Starting time of the MPI search run.

		:py:func:`time()`

	STDOUT : skysearch_lib.py
		
		Used to output progress messages and information to the terminal.

		:py:func:`os.sys.stdout.write`.
	.. <END - SKYSEARCH_LIB>

.. <DATA> 2018-03-24 #########################################################
.. glossary::

	Galaxia : Data

		The program responsible for producing the baryon data skysearcher is processing.

		`Galaxia Website <http://galaxia.sourceforge.net/>`_
	
	satid : Data

		A unique integer value representing the in falling satellite each star belongs to.

		This value comes from :term:`Galaxia`.

	m_book : Data

		A :py:mod:`numpy.ndarray` of satellite mass data.
		
	t_book : Data

		A :py:mod:`numpy.ndarray` of satellite accretion time data.
		
	j_book : Data

		A :py:mod:`numpy.ndarray` of satellite circularity data.

.. <SEARCH ALGORTHM> 2018-03-24 ##############################################
.. glossary::

	annulus : Search-algorithm

		A pair of inner and outer radial values which form a ring shaped area for searching within. 

	.. <HALO> ################################################################

	halo : Halo-specific-data

		The name of the halo.  "halo02"

		From :term:`Galaxia`
	
	local_satid_table : Halo-specific-data

		A sub-table of :term:`satid_table`. 

	satid_list : Halo-specific-data

		A :py:obj:`list` of all :term:`satid`'s belonging to the particular :term:`halo`.

		:func:`satid_setup`

	satid_table : Halo-specific-data

		An :py:mod:`astropy.table.Table` with columns for :term:`satid`, PHI and RADS.

	grids :  Halo-specific-data

		A :py:obj:`list` of tuples containing a halo name and an absolute PATH to the :py:mod:`numpy.ndarray`. 

		The result of :func:`grid_list`

		.. code-block:: python3

			[('halo07', '~/$PATH/skysearcher/data/grids/halo07_4.0Mpc_h158_grid.npy'),
			 ('halo09', '~/$PATH/skysearcher/data/grids/halo09_4.0Mpc_h158_grid.npy'),
			 ('halo20', '~/$PATH/skysearcher/data/grids/halo20_4.0Mpc_h158_grid.npy'),
			 ('halo05', '~/$PATH/skysearcher/data/grids/halo05_4.0Mpc_h158_grid.npy'),
			 ('halo17', '~/$PATH/skysearcher/data/grids/halo17_4.0Mpc_h158_grid.npy'),
			 ('halo10', '~/$PATH/skysearcher/data/grids/halo10_4.0Mpc_h158_grid.npy'),
			 ('halo14', '~/$PATH/skysearcher/data/grids/halo14_4.0Mpc_h158_grid.npy'),
			 ('halo12', '~/$PATH/skysearcher/data/grids/halo12_4.0Mpc_h158_grid.npy'),
			 ('halo08', '~/$PATH/skysearcher/data/grids/halo08_4.0Mpc_h158_grid.npy'),
			 ('halo15', '~/$PATH/skysearcher/data/grids/halo15_4.0Mpc_h158_grid.npy'),
			 ('halo02', '~/$PATH/skysearcher/data/grids/halo02_4.0Mpc_h158_grid.npy')]


	grid : Halo-specific-data
		
		A :py:mod:`numpy.ndarray` containing all the needed data regarding the projected stellar content of a :term:`halo`.

		The result of :func:`load_grid`

	.. <END HALO> ############################################################


	.. <FEATURE-PROPERTIES> 2018-03-24 #######################################
	
	r : Feature-property

		The center radial value for a feature in units of :term:`Kpc`.

	r_start : Feature-property

		The starting radius (always < :term:`r`) for a feature in units of :term:`Kpc`.

	r_stop : Feature-property

		The stopping radius (always > :term:`r`) for a feature in units of :term:`Kpc`.

	xbmax : Feature-property

		A :py:obj:`list` for the max values for a feature.

	log_n_stars_max : Feature-property

		A :py:obj:`list` for the max Log10(:term:`n_stars`) values for a feature.

	sat_stars : Feature-property

		A :py:obj:`dict` containing keys for :term:`satid` and values for number of stars.

		:func:`new_sat_stars`

	n_boxes : Feature-property

		The value for total array elements in a feature.

	n_stars : Feature-property

		The total number of stars in a feature.

	n_segments : Feature-property

		The run length (number of continuous regions) for a feature.

	starting_deg : Feature-property

		For a region, the starting angular value in :term:`Radians`.

	n_skips : Feature-property

		The allowed segment skips.  Should be set by :term:`N_SKIPS`.

		.. code-block:: python3
		
		    n_skips = N_SKIPS

	angular_extent : Feature-property

		The feature's total angular extent in :term:`Radians`.

	domsat_id : Feature-property

		An :py:obj:`int` value for the most abundant satellite id (:term:`satid`).

	.. <END FEATURE-PROPERTIES> ##############################################


	.. <SEARCH REGION> #######################################################
	
	_deg0 : Search-region

		The starting angle of a search region in :term:`Radians`.

	_deg1 : Search-region

		The final angle of a search region in :term:`Radians`.

	mu_2 : Search-region

		The mean number of stars in the local section of the current :term:`annulus`.

	r_idx2 : Search-region

		The array indices's for this sub-annulus-section which :term:`mu_2` is calculated.

	r_info : Search-region

		A :py:obj:`list` of the current search region's area information.

		.. code-block:: python3
		
		    r_info = [r_start, r_stop, _deg0, _deg1]


	.. <END SEARCH REGION> ###################################################
	
	n_stars_here : mpi_search.py

		Number of stars in this search region.

		.. code-block:: python3
		
			n_stars_here = grid[:, :, 1][idx].sum()

	xbox : mpi_search.py

		The array of contrast density values for each array element included in a given search region.

		:func:`get_xbox`

	n_boxes_in_seg : mpi_search.py

		The number of grid boxes (array elements) this search region covers.
	
	n_seg_increase : mpi_search.py

		Number of segments continuously increasing in the number of stars.

	n_seg_decrease : mpi_search.py

		Number of segments continuously decreasing in the number of stars.

	last_nstars : mpi_search.py

		The last count of number of stars.

	n_mt_seg : mpi_search.py
		
		Number of empty segments.

	a_tic : mpi_search.py

		Start time

	one_before : mpi_search.py

		Boolean value to indicate the existence of a immediately previous accepted segment.

	mod : mpi_search.py
		
		The result of :func:`kpc_to_arcmin` and serves as the ratio between :term:`Kpc` and :term:`Arc-min`.

	.. <END - SEARCH ALGORTHM>

	.. <SEARCH-EXTENT> #######################################################

	_radii : Search-extent
		
		:py:obj:`list` of radii to search from :term:`R_START` to :term:`R_STOP`.

		The result of :func:`radii`
	
	n_skips : Search-extent

		Integer value representing the number of remaining segment skips.

		:term:`N_SKIPS`

	.. <END SEARCH-EXTENT> ###################################################

	

	.. <FUNCTIONS>
	.. <END - FUNCTIONS>

:ref:`genindex`
