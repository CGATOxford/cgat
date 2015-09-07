.. _CGATInstallationLinux:

==================
Linux Installation
==================

Here we show different ways of installing the CGAT code on Linux
so you can take the one that better suits your requirements.

If you want to try out the CGAT code in a new, isolated 
development environment, we recommend you to go for the 
:ref:`AutomatedInstall`. 

In case you are using Conda already you may want to have a look
at the :ref:`CondaInstall` in order to integrate our code with yours.

Moreover, we also describe how to install the CGAT code with ``pip`` 
in section :ref:`PipInstall` but please note that this installation 
process is no longer maintained and may cause you problems due to
missing dependencies.

Introduction
============

All CGAT code is hosted at Git Hub::

        https://github.com/CGATOxford/cgat

The code contains a wide range of scripts and pipelines to analyse
Next-Generation Sequencing Data. The code in Git Hub is constantly
changing so it is considered a ``development`` version.

On the other hand, CGAT selected a smaller number of scripts (only 
scripts, no pipelines) that are being most actively used by the group. 
This is what we call the CGAT Code Collection and it was published 
in Bioinformatics_. The scripts included in the Code Collection can
be considered as a ``stable`` version and they are ready for production
out there. The CGAT Code Collection is available at PyPi::

        https://pypi.python.org/pypi/CGAT

CGAT code has many dependencies so installing it on different Linux
platforms from scratch is difficult by default. That is why we assist
the installation process with both Conda_ and an automated script.

Currently we provide four ways of getting the CGAT code installed, namely:

* ``cgat-scripts``

* ``cgat-scripts-lite``

* ``cgat-devel``

* ``cgat-devel-lite``

``cgat-scripts`` and ``cgat-scripts-lite`` include the CGAT Code 
Collection published in Bioinformatics_. ``cgat-devel`` and 
``cgat-devel-lite`` contains the dependencies to work with the
scripts and the pipelines hosted in Git Hub.

We use the term ``lite`` to denote a type of installation where part
of the system dependencies (C and C++ compilers, boost library, etc.) 
are delegated to the operating system package manager during the
installation process. You may have these dependencies already installed 
on your system (see :ref:`PipInstall` for a detailed list). Otherwise, 
you will need root permissions to install them. The automated script 
only helps you with the installation of these system packages on Scientific
Linux 6.x, CentOS 6.x and Ubuntu 12.x.

Optionally, you could go for a ``full`` Conda installation instead and 
this will ask Conda to install the system dependencies as well in
a new, fresh environment. You do not need root permissions to do this
and, more importantly, this installation is virtually portable to
any Linux (x86_64) distribution.

Next, we explain how to use the automated script. Afterwards we show the
steps to follow if you want to use ``conda`` directly. Finally we describe the
installation process with ``pip``.

.. _AutomatedInstall:

Automated Installation
======================

First of all please get a copy of the installation script::

        wget https://raw.github.com/CGATOxford/cgat/master/install-CGAT-tools.sh
        chmod +x install-CGAT-tools.sh

Install the lite version
------------------------

The ``lite`` installation is suitable when all the system dependencies listed
in :ref:`PipInstall` are already available on the target system. If not, you can
use the ``--install-os-packages`` option of the installation script as root (or sudo)::

        ./install-CGAT-tools.sh --install-os-packages

Then, type::

        ./install-CGAT-tools.sh --cgat-scripts --lite [--location /path/to/folder]

or::

        ./install-CGAT-tools.sh --cgat-devel --lite [--location /path/to/folder]

The code will be installed in ``$HOME/cgat-install`` by default but you can optionally
set another path with the ``--location`` option.

Install the full version
------------------------

The ``full`` installation should work on any Linux (x86_64) distribution even if
the system dependencies are not installed previously. Just type::

        ./install-CGAT-tools.sh --cgat-scripts --full [--location /path/to/folder]

or::

        ./install-CGAT-tools.sh --cgat-devel --full [--location /path/to/folder]

Again, use the ``--location`` option if you prefer a custom installation folder. 

Please note that ``full`` installation is the default one so if you omit ``--lite``
or ``--full`` the latter option will be included automatically.

Enable CGAT environment
-----------------------

After installation, you need to enable the CGAT environment as follows::

        source $LOCATION/conda-install/bin/activate <installation-option>

where ``$LOCATION`` is the installation folder specified above or ``$HOME/cgat-install``.
From now on, you just need to type::

        cgat --help

Disable CGAT environment
------------------------

CGAT environment can be disabled by doing::

        source deactivate

Update CGAT installation
------------------------

If you want to update a working installation of the CGAT code, please do::

        ./install-CGAT-tools.sh --update [--location /path/to/location]

Uninstall CGAT
--------------

All the software installed with the script is located under a single folder so
uninstalling the CGAT code involves removing it. However, you can also use the
script to uninstall it as follows::

        ./install-CGAT-tools.sh --uninstall [--location /path/to/location]

.. _CondaInstall:

Installation with Conda
=======================

This option is suitable when you are already using a Conda environment and 
you want to integrate the CGAT code with your other software. In that case, 
you only need to type::

        # add cgat channel
        conda config --add channels cgat

        # install CGAT code in Conda's root environment
        conda install <cgat-package>

where ``<cgat-package>`` can be one of these:

* ``cgat-scripts``

* ``cgat-scripts-lite``

* ``cgat-devel``

* ``cgat-devel-lite``

Please note that ``cgat-devel`` and ``cgat-devel-lite`` are Conda metapackages
and you will need to clone CGAT's Git Hub repository to actually get the code::

        git clone https://github.com/CGATOxford/cgat.git

.. _PipInstall:

Installation with ``pip``
==========================

Please note that this installation process is no longer supported
so you may experience problems due to missing dependencies.

Install RPM dependencies (Scientific Linux 6.x or CentOS 6.x)
-------------------------------------------------------------

First of all, make sure that you have the EPEL (Extra Packages for Enterprise Linux) repository installed::

        wget http://dl.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm
        rpm -i epel-release-6-8.noarch.rpm

Then, use ``yum`` to install all system dependencies::

        yum install gcc                 # required by python
        yum install zlib-devel          # required by virtualenv
        yum install openssl-devel       # required by pip
        yum install bzip2-devel         # required by bx-python
        yum install gcc-c++             # required by pybedtools
        yum install freetype-devel      # required by matplotlib
        yum install libpng-devel        # required by matplotlib
        yum install blas atlas lapack   # required by scipy
        yum install gcc-gfortran        # required by scipy
        yum install postgresql-devel    # required by PyGreSQL
        yum install R-core-devel        # required by rpy2
        yum install readline-devel      # required by rpy2
        yum install mysql-devel         # required by MySQL-python
        yum install boost-devel         # required by alignlib
        yum install sqlite-devel        # required by CGAT

Now type this additional commands to get scipy working::

        ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
        ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so # Scientific Linux Only
        ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so


Install DEB dependencies (Ubuntu 12.04 LTS)
-------------------------------------------

Use ``apt-get`` to install all system dependencies::

        apt-get install gcc                  # required by python
        apt-get install zlib1g-dev           # required by virtualenv
        apt-get install libssl-dev           # required by pip
        apt-get install libbz2-dev           # required by bx-python
        apt-get install g++                  # required by pybedtools
        apt-get install libfreetype6-dev     # required by matplotlib
        apt-get install libpng12-dev         # required by matplotlib
        apt-get install libblas-dev          # required by scipy
        apt-get install libatlas-dev         # required by scipy
        apt-get install liblapack-dev        # required by scipy
        apt-get install gfortran             # required by scipy
        apt-get install libpq-dev            # required by PyGreSQL
        apt-get install r-base-dev           # required by rpy2
        apt-get install libreadline-dev      # required by rpy2
        apt-get install libmysqlclient-dev   # required by MySQL-python
        apt-get install libboost-dev         # required by alignlib
        apt-get install libsqlite3-dev       # required by CGAT

Build Python 2.7
----------------

Download and build your own, isolated Python installation::

        cd
        mkdir CGAT
        wget http://www.python.org/ftp/python/2.7.6/Python-2.7.6.tgz
        tar xzvf Python-2.7.6.tgz
        rm Python-2.7.6.tgz
        cd Python-2.7.6
        ./configure --prefix=$HOME/CGAT/Python-2.7.6
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
.. _Bioinformatics: http://www.ncbi.nlm.nih.gov/pubmed/24395753
