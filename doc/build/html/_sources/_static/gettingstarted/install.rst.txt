:ref:`genindex`

Install skysearcher
-------------------
The first thing is to run either `pip <https://docs.python.org/3/installing/index.html>`_ or `setup.py <https://docs.python.org/2/distutils/setupscript.html>`_ from skysearcher's top level (*root*) directory.  

Using `pip <https://docs.python.org/3/installing/index.html>`_
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Move into skysearcher's top level (*root*) directory:

.. code-block:: bash

    $USER@machine$ cd skysearcher

From skysearcher top directory:

.. code-block:: bash

    $USER@machine$ ls
    data      doc      paper      __pycache__      setup.py      skysearcher

Now ask `pip <https://docs.python.org/3/installing/index.html>`_ to install skysearcher:

.. code-block:: bash

    $USER@machine$ pip install ./ -v
    Processing ~/$USER/skysearcher
    Requirement already satisfied: numpy in ~/$USER/lib/python3.6/site-packages (from skysearcher==0.0.1)
    Requirement already satisfied: mpi4py in ~/$USER/lib/python3.6/site-packages (from skysearcher==0.0.1)
    Requirement already satisfied: numba in ~/$USER/lib/python3.6/site-packages (from skysearcher==0.0.1)
    Requirement already satisfied: astropy in ~/$USER/lib/python3.6/site-packages (from skysearcher==0.0.1)
    Requirement already satisfied: matplotlib in ~/$USER/lib/python3.6/site-packages (from skysearcher==0.0.1)
    Requirement already satisfied: llvmlite in ~/$USER/lib/python3.6/site-packages (from numba->skysearcher==0.0.1)
    Requirement already satisfied: pytest>=2.8 in ~/$USER/lib/python3.6/site-packages (from astropy->skysearcher==0.0.1)
    Requirement already satisfied: six>=1.10 in ~/$USER/lib/python3.6/site-packages (from matplotlib->skysearcher==0.0.1)
    Requirement already satisfied: python-dateutil>=2.0 in ~/$USER/lib/python3.6/site-packages (from matplotlib->skysearcher==0.0.1)
    Requirement already satisfied: pytz in ~/$USER/lib/python3.6/site-packages (from matplotlib->skysearcher==0.0.1)
    Requirement already satisfied: cycler>=0.10 in ~/$USER/lib/python3.6/site-packages (from matplotlib->skysearcher==0.0.1)
    Requirement already satisfied: pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 in ~/$USER/lib/python3.6/site-packages (from matplotlib->skysearcher==0.0.1)
    Requirement already satisfied: py>=1.5.0 in ~/$USER/lib/python3.6/site-packages (from pytest>=2.8->astropy->skysearcher==0.0.1)
    Requirement already satisfied: setuptools in ~/$USER/lib/python3.6/site-packages (from pytest>=2.8->astropy->skysearcher==0.0.1)
    Requirement already satisfied: attrs>=17.2.0 in ~/$USER/lib/python3.6/site-packages (from pytest>=2.8->astropy->skysearcher==0.0.1)
    Requirement already satisfied: pluggy<0.7,>=0.5 in ~/$USER/lib/python3.6/site-packages (from pytest>=2.8->astropy->skysearcher==0.0.1)
    Installing collected packages: skysearcher
        Running setup.py install for skysearcher ... done
    Successfully installed skysearcher-0.0.1
    You are using pip version 9.0.1, however version 9.0.2 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.      

Using `setup.py <https://docs.python.org/2/distutils/setupscript.html>`_
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Move into skysearcher's top (*root*) directory:

.. code-block:: bash

    $USER@machine$ cd skysearcher

From skysearcher top directory:

.. code-block:: bash

    $USER@machine$ ls
    data      doc      paper      __pycache__      setup.py      skysearcher

Ask `setup.py <https://docs.python.org/2/distutils/setupscript.html>`_ to build skysearcher:

.. code-block:: bash

    $USER@machine$ python setup.py build
    running build
    running build_py
    creating build
    creating build/lib
    creating build/lib/skysearcher
    copying skysearcher/mpi_search.py -> build/lib/skysearcher
    copying skysearcher/__init__.py -> build/lib/skysearcher
    copying skysearcher/skysearch_lib.py -> build/lib/skysearcher
    copying skysearcher/new_cfg.py -> build/lib/skysearcher
    running egg_info
    creating skysearcher.egg-info
    writing skysearcher.egg-info/PKG-INFO
    writing dependency_links to skysearcher.egg-info/dependency_links.txt
    writing requirements to skysearcher.egg-info/requires.txt
    writing top-level names to skysearcher.egg-info/top_level.txt
    writing manifest file 'skysearcher.egg-info/SOURCES.txt'
    reading manifest file 'skysearcher.egg-info/SOURCES.txt'
    writing manifest file 'skysearcher.egg-info/SOURCES.txt'    

Ask `setup.py <https://docs.python.org/2/distutils/setupscript.html>`_ to install skysearcher:

.. code-block:: bash

    $USER@machine$ python setup.py install
    running install
    running bdist_egg
    running egg_info
    writing skysearcher.egg-info/PKG-INFO
    writing dependency_links to skysearcher.egg-info/dependency_links.txt
    writing requirements to skysearcher.egg-info/requires.txt
    writing top-level names to skysearcher.egg-info/top_level.txt
    reading manifest file 'skysearcher.egg-info/SOURCES.txt'
    writing manifest file 'skysearcher.egg-info/SOURCES.txt'
    installing library code to build/bdist.linux-x86_64/egg
    running install_lib
    running build_py
    creating build/bdist.linux-x86_64
    creating build/bdist.linux-x86_64/egg
    creating build/bdist.linux-x86_64/egg/skysearcher
    copying build/lib/skysearcher/mpi_search.py -> build/bdist.linux-x86_64/egg/skysearcher
    copying build/lib/skysearcher/__init__.py -> build/bdist.linux-x86_64/egg/skysearcher
    copying build/lib/skysearcher/skysearch_lib.py -> build/bdist.linux-x86_64/egg/skysearcher
    copying build/lib/skysearcher/new_cfg.py -> build/bdist.linux-x86_64/egg/skysearcher
    byte-compiling build/bdist.linux-x86_64/egg/skysearcher/mpi_search.py to mpi_search.cpython-36.pyc
    byte-compiling build/bdist.linux-x86_64/egg/skysearcher/__init__.py to __init__.cpython-36.pyc
    byte-compiling build/bdist.linux-x86_64/egg/skysearcher/skysearch_lib.py to skysearch_lib.cpython-36.pyc
    byte-compiling build/bdist.linux-x86_64/egg/skysearcher/new_cfg.py to new_cfg.cpython-36.pyc
    creating build/bdist.linux-x86_64/egg/EGG-INFO
    copying skysearcher.egg-info/PKG-INFO -> build/bdist.linux-x86_64/egg/EGG-INFO
    copying skysearcher.egg-info/SOURCES.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
    copying skysearcher.egg-info/dependency_links.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
    copying skysearcher.egg-info/requires.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
    copying skysearcher.egg-info/top_level.txt -> build/bdist.linux-x86_64/egg/EGG-INFO
    zip_safe flag not set; analyzing archive contents...
    creating dist
    creating 'dist/skysearcher-0.0.1-py3.6.egg' and adding 'build/bdist.linux-x86_64/egg' to it
    removing 'build/bdist.linux-x86_64/egg' (and everything under it)
    Processing skysearcher-0.0.1-py3.6.egg
    Copying skysearcher-0.0.1-py3.6.egg to ~/$USER/lib/python3.6/site-packages
    Adding skysearcher 0.0.1 to easy-install.pth file

    Installed ~/$USER/lib/python3.6/site-packages/skysearcher-0.0.1-py3.6.egg
    Processing dependencies for skysearcher==0.0.1
    Searching for matplotlib==2.1.1
    Best match: matplotlib 2.1.1
    Adding matplotlib 2.1.1 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for astropy==2.0.3
    Best match: astropy 2.0.3
    Adding astropy 2.0.3 to easy-install.pth file
    Installing fits2bitmap script to ~/$USER/bin
    Installing fitscheck script to ~/$USER/bin
    Installing fitsdiff script to ~/$USER/bin
    Installing fitsheader script to ~/$USER/bin
    Installing fitsinfo script to ~/$USER/bin
    Installing samp_hub script to ~/$USER/bin
    Installing volint script to ~/$USER/bin
    Installing wcslint script to ~/$USER/bin

    Using ~/$USER/lib/python3.6/site-packages
    Searching for numba==0.36.2
    Best match: numba 0.36.2
    Adding numba 0.36.2 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for mpi4py==2.0.0
    Best match: mpi4py 2.0.0
    Adding mpi4py 2.0.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for numpy==1.13.3
    Best match: numpy 1.13.3
    Adding numpy 1.13.3 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for pyparsing==2.2.0
    Best match: pyparsing 2.2.0
    Adding pyparsing 2.2.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for cycler==0.10.0
    Best match: cycler 0.10.0
    Adding cycler 0.10.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for pytz==2017.3
    Best match: pytz 2017.3
    Adding pytz 2017.3 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for python-dateutil==2.6.1
    Best match: python-dateutil 2.6.1
    Adding python-dateutil 2.6.1 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for six==1.11.0
    Best match: six 1.11.0
    Adding six 1.11.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for pytest==3.3.0
    Best match: pytest 3.3.0
    Adding pytest 3.3.0 to easy-install.pth file
    Installing py.test script to ~/$USER/bin
    Installing pytest script to ~/$USER/bin

    Using ~/$USER/lib/python3.6/site-packages
    Searching for llvmlite==0.21.0
    Best match: llvmlite 0.21.0
    Adding llvmlite 0.21.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for pluggy==0.6.0
    Best match: pluggy 0.6.0
    Adding pluggy 0.6.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for attrs==17.3.0
    Best match: attrs 17.3.0
    Adding attrs 17.3.0 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Searching for setuptools==36.5.0.post20170921
    Best match: setuptools 36.5.0.post20170921
    Adding setuptools 36.5.0.post20170921 to easy-install.pth file
    Installing easy_install script to ~/$USER/bin
    Installing easy_install-3.6 script to ~/$USER/bin

    Using ~/$USER/lib/python3.6/site-packages
    Searching for py==1.5.2
    Best match: py 1.5.2
    Adding py 1.5.2 to easy-install.pth file

    Using ~/$USER/lib/python3.6/site-packages
    Finished processing dependencies for skysearcher==0.0.1

:ref:`genindex`
