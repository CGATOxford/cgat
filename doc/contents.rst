.. Gpipe documentation master file, created by
   sphinx-quickstart on Mon Jul 20 11:17:20 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the CGAT script collection
=====================================

This document contains documentation for scripts and libraries written
for/by CGAT. 

.. note::
   The documentation is currently under construction.

.. toctree::
   :maxdepth: 2

   scripts.rst
   modules.rst

Installation
============

TODO

Disclaimer
==========

This collection of scripts is the result of 10 years working in various 
fields in bioinformatics. It contains both the good, the bad and the ugly. 
Use at your own risk.

Other projects
==============

There are a couple of separate libraries that
have python bindings and are used by the 
script collection.

These are:

   * alignlib_: a C++ library for sequence alignment with python bindings
   * ncl_: an implementation of nested containment lists
   * fastgtf_: a fast parser for :term:`GTF` formatted files
   * SphinxReport: a report generator based on Sphinx
   * pysam: python interface for SAM/BAM-files
   * Components: library for computing connected components quickly

Pipelines
=========

.. toctree::
   :maxdepth: 2

   pipelines.rst
   setup.rst
   BuildingPipelines.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _alignlib: http://wwwfgu.anat.ox.ac.uk/~andreas/alignlib
.. _ncl: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/ncl/contents.html
.. _fastgtf: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/fastgtf/contents.html
