.. image:: https://travis-ci.org/CGATOxford/cgat.svg?branch=master
    :target: https://travis-ci.org/CGATOxford/cgat

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io/recipes/cgat-scripts/README.html

===========================
The CGAT code collection
===========================

The CGAT Code collection has two components. The first component
is a collection of scripts, the CGAT tools. Scripts are located
`here <https://github.com/CGATOxford/cgat/tree/master/CGAT/scripts>`_
and can be ran on the commandline. Within this repository we also have a
number of utility functions that help with working with varuious file formats
in python. These are located `here <https://github.com/CGATOxford/cgat/tree/master/CGAT>`_.

The second component is a collection of pipelines that utilise the
functionality of the CGAT tool collection and these pipelines can be accessed
`here <https://github.com/CGATOxford/CGATPipelines>`_.

For questions, please subscribe and contact us at the 
`CGAT user group
<https://groups.google.com/forum/?fromgroups#!forum/cgat-user-group>`_.

Documentation of CGAT tools is available 
`here <https://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html>`_.


CGAT Tools
==========

Install using Conda
-------------------

The preferred method to install the CGAT code collection is using the installation script, which uses conda_.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/CGATOxford/cgat/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install set of production scripts (well tested):
        bash install-CGAT-tools.sh --production [--location </full/path/to/folder/without/trailing/slash>]

        # or go for the latest development version:
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

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

Example usage of a script
-------------------------

To use CGAT Tools, use the ``cgat`` front-end. For example, to
strip sequence and quality information from a bam_ file, type:

   cgat bam2bam --strip=sequence < in.bam > out.bam

For more extensive examples please refer to the documentation 
`here <https://www.cgat.org/downloads/public/cgat/documentation/CGATReference.html>`_


CGAT Pipelines
==============

We have developed numerous pipelines in comparative genomics
and NGS analysis. The pipelines are generally available and should
be fairly portable. Some documentation of the pipelines is 
`here <https://www.cgat.org/downloads/public/cgatpipelines/documentation/>`_.

They are under continuous development and changing rapidly.
However, they might give some ideas or building blocks when developing
your own pipelines.

.. _bam: http://en.wikipedia.org/wiki/SAMtools
.. _CGAT: http://www.cgat.org


