.. _CleanInstall:

============================
CGAT Code Clean Installation
============================

The section describes the steps to install the CGAT Code collection and its
dependencies inside a newly created environment. The instructions below have
been tested on a Red Hat Linux Enterprise 6.x-based operating system.
For instructions for OS X, see :ref:`CGATInstallationOSX`.

The :ref:`CleanQuickInstallation` uses CGAT supplied scripts for
installation, while :ref:`CleanManualInstallation` lists all the 
steps individually.

.. _CleanQuickInstallation:

Quick installation
==================

Get a copy of the installation scripts
--------------------------------------

Download and place them into your home directory::

        cd
        wget https://raw.github.com/CGATOxford/cgat/master/requires.txt
        wget https://raw.github.com/CGATOxford/cgat/master/setup-RPMs.sh
        wget https://raw.github.com/CGATOxford/cgat/master/setup-CGAT.sh

Install RPM dependencies
------------------------

Become root (or ask your system administrator to do it for you) and run ``setup-RPMs.sh``::

        ./setup-RPMs.sh

Install a Python virtual environment with the CGAT code collection
------------------------------------------------------------------- 

Do not be root for this step and run ``setup-CGAT.sh``::

        ./setup-CGAT.sh

Test the installation
---------------------

First activate the CGAT virtual environment::

        source $HOME/CGAT/virtualenv-1.10.1/cgat-venv/bin/activate

Then, test the cgat command::

        cgat --help

Finish the CGAT virtual environment
-----------------------------------

When you are done, you may deactivate the CGAT virtual environment::

        deactivate


.. _CleanManualInstallation:

Manual installation
===================

Install RPM dependencies
------------------------

You can either install them one by one or all at the same time with yum::

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
        yum install sqlite-devel        # requited by CGAT

and perform the additional configuration for scipy::

        ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
        ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so
        ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so

Build Python 2.7.5
------------------

Download and build your own, isolated Python installation::

        cd
        mkdir CGAT
        wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
        tar xzvf Python-2.7.5.tgz
        rm Python-2.7.5.tgz
        cd Python-2.7.5
        ./configure --prefix=$HOME/CGAT/Python-2.7.5
        make
        make install
        cd
        rm -rf Python-2.7.5

Create a virtual environment
----------------------------

Create an isolated virtual environment where all your Python packages will be installed::

        cd
        cd CGAT
        curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
        tar xvfz virtualenv-1.10.1.tar.gz
        rm virtualenv-1.10.1.tar.gz
        cd virtualenv-1.10.1
        $HOME/CGAT/Python-2.7.5/bin/python virtualenv.py cgat-venv
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


