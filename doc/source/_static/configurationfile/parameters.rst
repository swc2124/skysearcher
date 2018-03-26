:ref:`genindex`

Configuration parameters
------------------------
.. note:: All ``PATH`` parameters are absolute.

.. _config-PATH:

PATH
++++
``DATA_DIR`` : str
	Path to data directory
``PLOT_DIR`` : str
    Path to plot output directory.
``TABLE_DIR`` : str
    Path to hdf5 table input and output directory.
``MPI_TABLE_DIR`` : str
    Path to mpi hdf5 table directory within :term:`TABLE_DIR`.
``GRID_DIR`` : str
    Path to numpy array input files directory.
``TABLE_FORMAT`` : str
    Table format for astropy.tables.Table (default: {"hdf5"})
``TABLE_HDF5_PATH`` : str
    Table path for hdf5 table (default: {"hdf5"})
``TABLE_EXT`` : str
    Table file-type extension (default: {"hdf5"})
``GRID_EXT`` : str
    Extension for array files. (default: {"npy"})

.. _config-Search Extent:

Search Extent
+++++++++++++
``R_START`` : int
    Starting radius in units of :term:`Kpc` to start search from.
``R_STOP`` : int
    Ending search radius in units of :term:`Kpc`.
``R_STEP`` : int
    Step interval in units of :term:`Kpc`
``R_SCALE`` : float
    Percent to multiply radius by to get + and - :term:`Kpc` for annulus.
``ANNULUS_PHI_STEP`` : float
    How many sections to divide each annulus by.

.. _config-Accept Reject:

Accept Reject
+++++++++++++
``XBOX_CUT`` : float
    Minimum contrast density a feature must have.
``MIN_LOG_NSTARS`` : float
    Minimum log10 value for number of stars a feature must have.
``MIN_N_SEGMENTS`` : int
    Minimum number of segments a feature must be
``N_SKIPS`` : int
    Allowed number of skipped segments before ending feature.

.. _config-Run Time:

Run Time
++++++++
``SAVE_INTERVAL`` : int
    Number of features to hold between saves.

.. _config-DATA:

DATA
++++
``D_MPC`` : float
    Distance to target halo in units of :term:`Mpc`.

:ref:`genindex`
