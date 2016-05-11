.. _CGATInstallationSL6:

=================================
Scientific Linux 6.X Installation
=================================

Scientific Linux is a Red Hat Enterprise Linux-based operating 
system and therefore the instructions below should apply to all
other derivatives as well.

The :ref:`SL6QuickInstallation` uses CGAT supplied scripts for
installation, while :ref:`SL6ManualInstallation` lists all the 
steps individually.

.. _SL6QuickInstallation:

Quick installation
==================

Get a copy of the installation scripts
--------------------------------------

Download and place them into your home directory::

        cd
        wget https://raw.github.com/CGATOxford/cgat/master/install-CGAT-tools.sh
        chmod +x install-CGAT-tools.sh

Install RPM dependencies
------------------------

Become root (or ask your system administrator to do it for you) and run::

        ./install-CGAT-tools.sh --install-os-packages

Install the CGAT Code Collection
--------------------------------

CGAT uses Conda_ to create a new environment with all necessary software
pre-installed::

        ./install-CGAT-tools.sh --cgat-scripts

You do not need to be root for this step and the location of the installation
can be chosen with the ``--location <path>`` option.

Install a Python virtual environment with the CGAT code collection
------------------------------------------------------------------- 

Do not be root for this step and run::

        ./install-CGAT-tools.sh --install-python-deps

Test the installation
---------------------

First activate the CGAT virtual environment::

        source $HOME/CGAT-DEPS/virtualenv-1.11.6/cgat-venv/bin/activate

Then, test the cgat command::

        cgat --help

Finish the CGAT virtual environment
-----------------------------------

When you are done, you may deactivate the CGAT virtual environment::

        deactivate


.. _SL6ManualInstallation:

Manual installation
===================

Install RPM dependencies
------------------------

You can either install them one by one or all at the same time with ``yum``::

        yum install gcc                 # required by python
        yum install zlib-devel          # required by virtualenv
        yum install openssl-devel       # required by pip
        yum install bzip2-devel         # required by bx-python
        yum install gcc-c++             # required by pybedtools
        yum install freetype-devel      # required by matplotlib
        yum install libpng-devel        # required by matplotlib
        yum install blas atlas lapack   # required by scipy
        yum install gcc-gfortran        # required by scipyi
        yum install postgresql-devel    # required by PyGreSQL
        yum install R-core-devel        # required by rpy2
        yum install readline-devel      # required by rpy2
        yum install mysql-devel         # required by MySQL-python
        yum install boost-devel         # required by alignlib
        yum install sqlite-devel        # required by CGAT
        yum install mercurial           # required by bx-python

Please note that you may also need the EPEL (Extra Packages for Enterprise Linux) repository to install R::

        yum install epel-release

Now type this additional commands to get scipy working::

        ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
        ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so
        ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so

Build Python 2.7.6
------------------

Download and build your own, isolated Python installation::

        cd
        mkdir CGAT
        wget http://www.python.org/ftp/python/2.7.6/Python-2.7.6.tgz
        tar xzvf Python-2.7.6.tgz
        rm Python-2.7.6.tgz
        cd Python-2.7.6
        ./configure --column-prefix=$HOME/CGAT/Python-2.7.6
        make
        make install
        cd
        rm -rf Python-2.7.6

Create a virtual environment
----------------------------

Create an isolated virtual environment where all your Python packages will be installed::

        cd
        cd CGAT
        wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz
        tar xvfz virtualenv-1.11.6.tar.gz
        rm virtualenv-1.11.6.tar.gz
        cd virtualenv-1.11.6
        $HOME/CGAT/Python-2.7.6/bin/python virtualenv.py cgat-venv
        source cgat-venv/bin/activate

Install Python dependencies
---------------------------

Use pip to install all the packages on which CGAT Code Collection depends on::

        pip install cython
        pip install numpy
        pip install pysam
        pip install https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2
        pip install biopython
        pip install pybedtools
        pip install matplotlib
        pip install scipy
        pip install -r https://raw.github.com/CGATOxford/cgat/master/requires.txt
        pip install CGAT

Test CGAT Code Collection
-------------------------

If everything went fine with the previous steps you should be able to execute
the following command::

        cgat --help

Finish the CGAT virtual environment
-----------------------------------

When you are done, you may deactivate the CGAT virtual environment::

        deactivate

.. _Conda: http://conda.pydata.org/
