.. image:: https://travis-ci.org/CGATOxford/cgat.svg?branch=master
    :target: https://travis-ci.org/CGATOxford/cgat

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io/recipes/cgat-scripts/README.html

===========================
The CGAT Code Collection
===========================

The CGAT_ Code Collection has two components. The first component
is a collection of scripts in this repository, which are located
`here <https://github.com/CGATOxford/cgat/tree/master/CGAT/scripts>`_
and can be run using the ``cgat`` command. Within this repository we also have a
number of utility modules that help working with various file formats
in Python. These are located `here <https://github.com/CGATOxford/cgat/tree/master/CGAT>`_.

The second component is a collection of pipelines that utilise the
functionality of the scripts and can be accessed
`here <https://github.com/CGATOxford/CGATPipelines>`_.

For questions, please open a discussion on the GitHub 
`issue page <https://github.com/CGATOxford/cgat/issues>`_.

Documentation of CGAT tools is available 
`here <https://www.cgat.org/downloads/public/cgat/documentation/>`_.

Installation
============

Install using Conda
-------------------

The preferred method to install the CGAT code collection is using the installation script, which uses
`Conda <https://conda.io>`_.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/CGATOxford/cgat/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install set of production scripts (well tested):
        bash install-CGAT-tools.sh --production [--location </full/path/to/folder/without/trailing/slash>]

        # or go for the latest development version:
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # enable the conda environment as requested by the installation script:
        source </full/path/to/folder/without/trailing/slash>/conda-install/bin/activate cgat-s

        # finally, please run the cgatflow command-line tool to check the installation:
        cgat --help

The installation script will put everything under the specified location. The aim of the script is to
provide a portable installation that does not interfere with the existing software. As a result, you
will have a conda environment working with the CGAT scripts which can be enabled on demand according 
to your needs.

Install using pip
-----------------

You can also use pip to install the CGAT scripts. To go down this route, please type::

   pip install cgat

However, CGAT depends on numerous other python packages which themselves might require
manual intervention. Therefore, our preferred method of installation is through conda. 

Usage
=====

Run the ``cgat --help`` command to see what scripts are available and how to use them.
For example, to strip sequence and quality information from a bam_ file, type::

   cgat bam2bam --strip=sequence < in.bam > out.bam

For more extensive examples please refer to the documentation 
`here <https://www.cgat.org/downloads/public/cgat/documentation/CGATReference.html>`_

.. _bam: http://en.wikipedia.org/wiki/SAMtools
.. _CGAT: http://www.cgat.org
