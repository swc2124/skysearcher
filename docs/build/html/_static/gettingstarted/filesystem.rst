:ref:`genindex`

The file system
------------------------------

.. parsed-literal::

    **skysearcher/**                       <-- root directory
        │
        ├── ``setup.py``                    <-- setup.py script
        │
        ├── **data/**                      <-- directory for all data files
        │       │
        │       ├── **grids/**             <-- numpy arrays used by mpi_search.py
        │       ├── **plots/**             <-- output plots
        │       ├── **tables/**            <-- temp hdf5 tables used by MPI
        │       │
        │       ├── satj_array.npy     <-- satellite circularity
        │       ├── satmass_array.npy  <-- satellite mass
        │       ├── satprop.ebf        <-- satellite properties
        │       └── satage_array.npy   <-- satellite age
        │
        └── **skysearcher/**               <-- package
                │
                ├── __init__.py        <-- init.py
                ├── **rc.cfg**             <-- *main config file*
                ├── new_cfg.py         <-- for making new config file
                ├── mpi_search.py      <-- main loop
                └── skysearch_lib.py   <-- all functions

:ref:`genindex`
