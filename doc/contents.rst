.. Gpipe documentation master file, created by
   sphinx-quickstart on Mon Jul 20 11:17:20 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================================
Welcome to the CGAT code collection
=====================================

This document brings together the various pipelines and scripts
written before and during CGAT.

.. note::
   The documentation is under construction.

Overview
========

The CGAT code collection has grown out of the work in comparative genomics
by the Ponting group in the last decade. Now, CGAT has added
functionality to do next-generation sequencing analysis. 

The collection has three major components, these are directories in 
the package.

* scripts
  A collection of useful scripts for genomics and NGS analysis

* CGAT
  A collection of modules with utility functions for genomics and NGS analysis.
  
* CGATPipelines
  A collection of pipelines for common workflows in genomics and NGS analysis.

Pipelines
=========

CGAT has written a few :term:`production pipelines`. These pipelines perform
basic tasks, are fairly generic and might be of wider interest. Below is a 
complete list:
 
.. toctree::
   :maxdepth: 1	

   scripts/pipeline_annotations.rst   
   scripts/pipeline_ancestral_repeats.rst
   scripts/pipeline_chains.rst
   scripts/pipeline_chipseq.rst
   scripts/pipeline_liftover.rst
   scripts/pipeline_readqc.rst
   scripts/pipeline_rnaseq.rst
   scripts/pipeline_variants.rst
   scripts/pipeline_benchmark_rnaseqmappers.rst
   scripts/pipeline_testing.rst
   scripts/pipeline_mapping.rst
   scripts/pipeline_rnaseqtranscripts.rst
   scripts/pipeline_transcriptome.rst
   scripts/pipeline_rnaseqdiffexpression.rst
   scripts/pipeline_peakcalling.rst
   scripts/pipeline_rnaseqlncrna.rst   

Help on installing, using and building pipelines is below:

.. toctree::
   :maxdepth: 2

   setup.rst
   pipelines.rst
   BuildingPipelines.rst
   PipelinesBackground.rst   
   MakePipelines.rst   

Scripts and modules
===================

The pipelines make use of a large number of scripts, tools and modules. Most
of these are written in python. The tools are grouped by topic:

.. toctree::
   :maxdepth: 2

   scripts.rst
   modules.rst

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

.. _alignlib: http://wwwfgu.anat.ox.ac.uk/~andreas/alignlib
.. _ncl: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/ncl/contents.html
.. _fastgtf: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/fastgtf/contents.html
