.. _contents:

=====================================
Welcome to the CGAT code collection
=====================================

This document brings together the various pipelines and scripts
written before and during CGAT. This documentation has two
parts. Below is the general documentation covering the complete
code collection. A different view on the code collection centering
on the scripts published in 
`Bioinformatics <http://www.ncbi.nlm.nih.gov/pubmed/24395753>`_ is
:ref:`here <cgat>`.

Overview
========

The CGAT code collection has grown out of the work in comparative genomics
by the Ponting group in the last decade. Now, CGAT has added
functionality to do next-generation sequencing analysis. 

The collection has three major components, these are directories in 
the package.

* :ref:`scripts`
  A collection of useful scripts for genomics and NGS analysis

* :ref:`modules`
  A collection of modules with utility functions for genomics and NGS analysis.
  
* :ref:`pipelines`
  A collection of pipelines for common workflows in genomics and NGS analysis.

Scripts and modules
===================

The CGAT code collection is as set of tools and modules for genomics.
Most of these scripts are written in python. The tools are grouped by topic:

.. toctree::
   :maxdepth: 2

   scripts.rst
   modules.rst

.. _cgatpipelines:

CGAT Pipelines
==============

CGAT pipelines perform basic tasks, are fairly generic and might be of wider interest. 

.. toctree::
   :maxdepth: 2

   Pipelines.rst

Developer's guide
=================

.. toctree::
   :maxdepth: 2

   developing.rst
   testing.rst
   styleguide.rst
   documenting.rst
   pipelines.rst
   release.rst
   galaxy.rst

Glossary
========

.. toctree::
   :maxdepth: 2

   glossary.rst   

Disclaimer
==========

The collection of scripts and tools is the outcome of 10 years working in various 
fields in bioinformatics. It contains both the good, the bad and the ugly. 
Use at your own risk.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

