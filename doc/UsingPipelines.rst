====================
Using CGAT pipelines
====================

This section provides a tutorial-like introduction to CGAT pipelines.

Introduction
=============

A pipeline takes input data and performs a series of automated steps
(:term:`task`) on it to produce some output data.

Each pipeline is usually coupled with a :term:`CGATReport` document
to summarize and visualize the results.

It really helps if you are familiar with following:

   * the unix command line to run and debug the pipeline
   * python_ in order to understand what happens in the pipeline
   * ruffus_ in order to understand the pipeline code
   * sge_ in order to monitor your jobs
   * git_ in order to up-to-date code

.. _PipelineSettingUp:

Setting up a pipeline
======================

*Step 1*: Get the latest clone of the cgat script repository.  Check
that your computing environment is appropriate (see :ref:`CGATSetup`).
The directory in which the CGAT code repository is located is the
:term:`source directory`. It will be abbreviated ``<cgat>`` in the
following commands. The source directory will contain the pipeline
master script named :file:`pipeline_<name>.py`, the default
configuration files and all the helper scripts and libraries to run
the pipeline.

*Step 2*: Create a :term:`working directory` and enter it. For example::

   mkdir version1
   cd version1

The pipeline will live there and all subsequent steps should be executed 
from within this directory.

*Step 3*: Obtain and edit an initial configuration file. Ruffus pipelines are controlled
by a configuration file. A configuration file with all the default values can be 
obtained by running::

      python <cgat>/CGATPipelines/pipeline_<name>.py config

This will create a new :file:`pipeline.ini` file. **YOU MUST EDIT THIS FILE**.
The default values are likely to use the wrong genome or point to non-existing
locations of indices and databases. The configuration file should be well documented
and the format is simple. The documenation for the
`ConfigParser <http://docs.python.org/library/configparser.html>`_ python module 
contains the full specification.

*Step 4*: Add the input files. The required input is specific for each pipeline; read
the pipeline documentation to find out exactly which files are needed
and where they should be put. Commonly, a pipeline
works from input files linked into the :term:`working directory` and named
following pipeline specific conventions.

*Step 5*: You can check if all the external dependencies to tools and
R packages are satisfied by running::

      python <cgat>/CGATPipelines/pipeline_<name>.py check

See :ref:`ExternalDependencies` for more information.

.. _PipelineRunning:

Running a pipeline
===================

Pipelines are controlled by a single python script called :file:`pipeline_<name>.py`
that lives in the :term:`source directory`. Command line usage information is available
by running::

   python <cgat>/pipeline_<name>.py --help

The basic syntax for ``pipeline_<name>.py`` is::

   python <cgat>/pipeline_<name>.py [options] _COMMAND_

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

clone <srcdir>
   clone a pipeline from :file:`srcdir` into the current
   directory. Cloning attempts to conserve disk space by linking.

In case you are running a long pipeline, make sure you start it appropriately, for example::

   nice -19 nohup <cgat>/pipeline_<name>.py make full

This will keep the pipeline running if you close the terminal.

Troubleshooting
---------------

Many things can go wrong while running the pipeline. Look out for

   * bad input format. The pipeline does not perform sanity checks on the input format.
       If the input is bad, you might see wrong or missing results or an error message.
   * pipeline disruptions. Problems with the cluster, the file system or the controlling terminal 
       might all cause the pipeline to abort.
   * bugs. The pipeline makes many implicit assumptions about the input files and the programs it
       runs. If program versions change or inputs change, the pipeline might not be able to deal with it.
       The result will be wrong or missing results or an error message.

If the pipeline aborts, locate the step that caused the error by reading the logfiles and
the error messages on stderr (:file:`nohup.out`). See if you can understand the error and guess
the likely problem (new program versions, badly formatted input, ...). If you are able to fix 
the error, remove the output files of the step in which the error occured and restart the 
pipeline. Processing should resume at the appropriate point.

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

   git pull

The first command retrieves the latest changes from the master repository
and the second command updates your local version with these changes.

.. _PipelineReporting:

Building pipeline reports
================================

Some of the pipelines are associated with an automated report generator to display
summary information as a set of nicely formatted html pages. In order to
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

