=========
Pipelines
=========

The pipelines below are ready-made pipelines for certain simple tasks. Pipelines are controlled by a ``ruffus``
script, though there are some legacy pipelines written in ``make``.

Introduction
=============

A pipeline takes input data and performs a series of automated steps (:term:`task`) on it to 
produce some output data. 

Each pipeline is usually coupled with a :term:`SphinxReport` document to summarize and 
visualize the results.

It really helps if you are familiar with following:

   * the unix command line to run and debug the pipeline
   * python_ in order to understand what happens in the pipeline
   * ruffus_ in order to understand the pipeline code
   * sge_ in order to monitor your jobs
   * mercurial_ in order to up-to-date code

.. _python: http:www.python.org
.. _ruffus: http://code.google.com/p/ruffus/
.. _sge: http://wikis.sun.com/display/GridEngine/Home
.. _mercurial: http://mercurial.selenic.com/wiki/
.. _sphinxreport: http://code.google.com/p/sphinx-report/

.. _PipelineSettingUp:

Setting up a pipeline
======================

Before starting, check that your computing environment is appropriate
(see :ref:`CGATSetup`). Once all components are in place, setting up a 
pipeline involves the following steps:

*Step 1*: Get the latest clone of the cgat script repository::

   hg clone http://www.cgat.org/hg/cgat/ src

.. note:: 
   You need to have mercurial installed.

The directory :file:`src` is the :term:`source directory`. It will be abbreviated
``<src>`` in the following commands. This directory will contain the pipeline
master script named :file:`pipeline_<name>.py`, the default configuration files
and all the helper scripts and libraries to run the pipeline.

*Step 2*: Create a :term:`working directory` and enter it. For example::

   mkdir version1
   cd version1

The pipeline will live there and all subsequent steps should be executed 
from within this directory.

*Step 3*: Obtain and edit an initial configuration file. Ruffus pipelines are controlled
by a configuration file. A configuration file with all the default values can be 
obtained by running::

      python <src>/pipeline_<name>.py config

This will create a new :file:`pipeline.ini` file. **YOU MUST EDIT THIS FILE**.
The default values are likely to use the wrong genome or point to non-existing
locations of indices and databases. The configuration file should be well documented
and the format is simple. The documenation for the
`ConfigParser <http://docs.python.org/library/configparser.html>`_ python module 
contains the full specification.

*Step 4*: Add the input files. The required input is specific for each pipeline; read
the pipeline documentation to find out exactly which files are needed. Commonly, a pipeline
works from input files copied or linked into the :term:`working directory` and named
following pipeline specific conventions.

.. _PipelineRunning:

Running a pipeline
===================

Pipelines are controlled by a single python script called :file:`pipeline_<name>.py`
that lives in the :term:`source directory`. Command line usage information is available
by running::

   python <src>/pipeline_<name>.py --help

The basic syntax for :term:`pipeline_<name>.py` is::

   python <src>/pipeline_<name>.py [options] _COMMAND_

``COMMAND`` can be one of the following:

make <task>
   run all tasks required to build :term:`task`

show <task>
   show tasks required to build :term:`task` without executing them

plot <task>
   plot image (requires `inkscape <http://inkscape.org/>`_) of pipeline state for :term:`task`

touch <task>
   touch files without running :term:`task` or its pre-requisites. This sets the 
   timestamps for files in :term:`task` and its pre-requisites such that they will 
   seem up-to-date to the pipeline.

config
   write a new configuration file :file:`pipeline.ini` with default values. An existing 
   configuration file will not be overwritten.

In case you are running a long pipeline, make sure you start it appropriately, for example::

   nice -19 nohup <src>/pipeline_<name>.py make full

This will keep the pipeline running if you close the terminal.

Troubleshooting
---------------

Many things can go wrong while running the pipeline. Look out for

   * bad input format. The pipeline does not perform sanity checks on the input format.
       If the input is bad, you might see wrong or missing results or an error message.
   * pipeline disrutions. Problems with the cluster, the file system or the controlling terminal 
       might all cause the pipeline to abort.
   * bugs. The pipeline makes many implicit assumptions about the input files and the programs it
       runs. If program versions change or inputs change, the pipeline might not be able to deal with it.
       The result will be wrong or missing results or an error message.

If the pipeline aborts, locate the step that caused the error by reading the logfiles and
the error messages on stderr (:file:`nohup.out`). See if you can understand the error and guess
the likely problem (new program versions, badly formatted input, ...). If you are able to fix 
the error, remove the output files of the step in which the error occured and restart the 
pipeline. It should continue from the appropriate location.

.. note::
   Look out for upstream errors. For example, the pipeline might build a geneset filtering
   by a certain set of contigs. If the contig names do not match, the geneset will be empty,
   but the geneset building step might conclude successfully. However, you might get an error
   in any of the downstream steps complaining that the gene set is empty. To fix this, fix
   the error and delete the files created by the geneset building step and not just the step
   that threw the error.

Updating to the latest code version
-----------------------------------

To get the latest bugfixes, go into the :term:`source directory` and type::

   hg pull
   hg update

The first command retrieves the latest changes from the master repository
and the second command updates your local version with these changes.

.. _PipelineDocumentation:

Building pipeline documentation
================================

Some of the pipelines are associated with an automated report generator to display
summary information as a set of nicely formatted :term:`html` pages. In order to
build the documentation, drop the appropriate :file:`conf.py` and :file:`sphinxreport.ini`
configuration files into the :term:`working directory` and run the pipeline command:

   nice -19 pipeline_<name>.py make build_report

This will create the report from scratch in the current directory. The report can
be viewed opening the file :file:`<work>/report/html/contents.html` in your browser.

Sphinxreport is quite powerful, but also runs quite slowly on large projects that
need to generate a multitude of plots and tables. In order to speed up this process,
there are some advanced features that Sphinxreport offers:

   * caching of results
   * multiprocessing
   * incremental builds
   * separate build directory

Please see the sphinxreport_ documentation for more information.

Ruffus-based pipelines
========================

Ruffus based pipelines start with the :file:`pipeline_` prefix. Each pipeline consists
of two files, the actual pipeline (e.g., :file:`pipeline_rnaseq.py`) and an associated
configuration file with default values for pipeline parameters (e.g., :file:`pipeline_rnaseq.ini`).

.. toctree::
   :maxdepth: 1	

   scripts/pipeline_annotations.rst   
   scripts/pipeline_ancestral_repeats.rst
   scripts/pipeline_chipseq.rst
   scripts/pipeline_rnaseq.rst
   scripts/pipeline_transcriptome.rst

Make-based pipelines
======================

.. toctree::
   :maxdepth: 1

   pipelines/CompareTranscripts.rst
   pipelines/Gpipe.rst
   pipelines/MapTranscripts454.rst

Background
============

There really are two types of pipelines. In ``production pipelines`` the inputs are usually
the same every time the pipeline is run and the output is known beforehand. For example, 
read mapping and quality control is a typical pipeline. These pipelines can be well optimized
and can be re-used with little change in configuration.

``analysis pipelines`` control scientific analyses and are much more in a state of flux. 
Here, the input might change over time as the analysis expands and the output will change
with every new insight or new direction a project takes. It will be still a pipeline as long as
the output can be generated from the input without manual intervention. These pipelines leave
less scope for optimization compared to ``production pipelines`` and adapting a pipeline to
a new project will involve significant refactoring.

In CGAT, we are primarily concerned with ``analysis pipelines``, though we have some 
``production pipelines`` for common tasks.

There are several ways to build pipelines. For example, there are generic workflow
systems like `taverna <http://www.taverna.org.uk>`_ which even provide GUIs for connecting
tasks. A developer writes some glue code permitting the output of one application to
be used as input for another application. Also, there are specialized workflow systems 
for genomics, for example `galaxy <http://galaxy.psu.edu>`_, which allows you to save and share
analyses. New tools can be added to the system and new data imported easily for example
from the UCSC genome browser.

Flexibility
   There always new tools and insights. A pipeline should be ultimately 
   flexible and not constraining us in the things we can do.

Scriptability
   The pipeline should be scriptable, i.e, the whole pipeline can be run within
   another pipeline. Similarly, parts of a pipeline can be duplicated to process 
   several data streams in parallel. This is a crucial feature in genome studies
   as a single analysis will not permit making inferences by itself. For example,
   consider you find in ChIP-Seq data from a particular transcription factor that
   it binds frequently in introns. You will need to run the same analysis on 
   data from other transcription factors in order to assess if intronic binding is
   remarkable.

Reproducibility
   The pipeline is fully automated. The same inputs and configuration will produce
   the same outputs.

Reusability
   The pipeline should be able to be re-used on similar data, maybe only requiring 
   changes to a configuration file.

Archivability
   Once finished, the whole project should be able to archived without too many
   major dependencies on external data. This should be a simple process and hence
   all project data should be self-contained. It should not involve going through 
   various directories or databases to figure out which files and tables belong
   to a project or a project depends on.

There probably is not one toolset to satisfy all these criteria.. We use the following 
tools to build a pipeline:

   * ruffus_ to control the main computational steps
   * sqlite_ to store the results of the computational steps
   * sphinxreport_ to visualize the data in the sqlite database

