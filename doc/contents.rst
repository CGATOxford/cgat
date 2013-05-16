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

.. _cgatpipelines:
 
.. toctree::
   :maxdepth: 1	

   pipelines/pipeline_annotations.rst   
   pipelines/pipeline_ancestral_repeats.rst
   pipelines/pipeline_chains.rst
   pipelines/pipeline_chipseq.rst
   pipelines/pipeline_liftover.rst
   pipelines/pipeline_readqc.rst
   pipelines/pipeline_rnaseq.rst
   pipelines/pipeline_variants.rst
   pipelines/pipeline_benchmark_rnaseqmappers.rst
   pipelines/pipeline_testing.rst
   pipelines/pipeline_mapping.rst
   pipelines/pipeline_rnaseqtranscripts.rst
   pipelines/pipeline_transcriptome.rst
   pipelines/pipeline_rnaseqdiffexpression.rst
   pipelines/pipeline_peakcalling.rst
   pipelines/pipeline_rnaseqlncrna.rst   
   pipelines/pipeline_capseq.rst
   pipelines/pipeline_cufflinks_optimization.rst
   pipelines/pipeline_exome.rst
   pipelines/pipeline_expression.rst
   pipelines/pipeline_fastqToBigWig.rst
   pipelines/pipeline_fusion.rst
   pipelines/pipeline_intervals.rst
   pipelines/pipeline_mappability.rst
   pipelines/pipeline_mapping_benchmark.rst
   pipelines/pipeline_medip.rst
   pipelines/pipeline_polyphen.rst
   pipelines/pipeline_promotors.rst
   pipelines/pipeline_quickstart.rst
   pipelines/pipeline_transfacmatch.rst
   pipelines/pipeline_variant_annotation.rst

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

Deverlopers guide

.. toctree::
   :maxdepth: 2

   testing.rst
   styleguide.rst
   release.rst
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

.. _alignlib: http://wwwfgu.anat.ox.ac.uk/~andreas/alignlib
.. _ncl: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/ncl/contents.html
.. _fastgtf: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/fastgtf/contents.html
