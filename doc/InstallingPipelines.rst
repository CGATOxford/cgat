.. _CGATSetup:

=========================
Installing CGAT pipelines
=========================

The CGAT pipelines, scripts and libraries make several assumptions about
the computing environment. This section describes how to install the code
and set up your computing environment.

Downloading and installing the source code
==========================================

To obtain the latest code, check it out from the public git_
repository and activate it::

   git clone https://github.com/CGATOxford/cgat.git
   cd cgat
   python setup.py develop

Please see the installation instructions for the `CGAT Toolkit
<http://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html>`_
if you run into problems.

Once checked-out, you can get the latest changes via pulling::

   git pull 

Some scripts contain cython_ code that needs to be recompiled if the
script or the pysam_ installation has changed. To rebuild all scripts,
for example after updating the repository, type::

   python cgat/scripts/cgat_rebuild_extensions.py

Recompilation requires a C compiler to be installed. 

Setting up the computing environment
====================================

The pipelines assume that Sun Grid Engine has been installed. Other
queueing systems might work, but expect to be disappointed. The
pipeline is started on a :term:`submit host` assuming a default queue
``all.q``. Other queues can be specified on the command line, for
example::

    python cgat/CGATPipelines/pipeline_<name>.py --cluster-queue=medium_jobs.q

A pipeline might start up to ``-p/--multiprocess`` processes. Preferentially,
tasks are sent to the cluster, but for some tasks this is not possible. 
These might thus run on the :term:`submit host`, so make sure it is fairly powerful.

Pipelines expects that the :term:`working directory` is accessible with
the same path both from the submit and the :term:`execution host`. 

Software requirements
=====================

On top of pipeline specific bioinformatics software, CGAT pipelines
make use a variety of software. Unfortunately we can't support many
versions. The following table gives a list software that are currently
being required, installed and the location where the dependency arises:

.. How to create this table:
.. python scripts/cgat_list_dependencies.py | tab2rst

+-------------+------------+---------------+-----------+----------------------------------------+
|tool         |required    |installed      |is_required|locations                               |
+-------------+------------+---------------+-----------+----------------------------------------+
|BroadPeak    |>=1.0       |?              |False      |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|DESeq        |>=1.17      |1.17.0         |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|DESeq2       |>=1.5.62    |1.5.62         |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|GATK         |>=2.7       |2.7-2          |True       |CGATPipelines/pipeline_exome.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|Gemini       |>=?         |-              |True       |CGATPipelines/pipeline_exome.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|HiddenMarkov |>=1.8.0     |1.8-0          |True       |CGATPipelines/PipelineRnaseq.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|MASS         |>=7.3.34    |7.3-34         |True       |CGATPipelines/PipelineRnaseq.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|MEDIPS       |>=1.15.0    |1.15.0         |True       |CGATPipelines/PipelineWindows.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|RColorBrewer |>=1.0.5     |1.0-5          |True       |CGAT/Expression.py                      |
|             |            |               |           |CGATPipelines/PipelineRnaseq.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|SICER        |>=1.1       |?              |False      |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|bedtools     |>=2.21.0    |2.21.0         |True       |CGATPipelines/PipelinePeakcalling.py    |
|             |            |               |           |CGATPipelines/PipelineWindows.py        |
|             |            |               |           |CGATPipelines/pipeline_windows.py       |
|             |            |               |           |scripts/runSPP.py                       |
+-------------+------------+---------------+-----------+----------------------------------------+
|bismark      |>=0.12.5    |0.12.5         |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|bowtie       |>=1.0.0     |1.0.0          |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|bowtie2      |>=2.2.3     |2.2.3          |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|bwa          |>=0.7.8     |0.7.8-r455     |True       |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|cufflinks    |>=2.2.1     |2.2.1          |True       |CGATPipelines/PipelineMapping.py        |
|             |            |               |           |CGATPipelines/PipelineRnaseq.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|edgeR        |>=3.7.16    |3.7.16         |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|fastq-dump   |>=2.1.7     |2.1.7          |True       |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|fastqc       |>=0.9.2     |0.9.2          |True       |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|featureCounts|>=1.4.3     |1.4.3-p1       |True       |CGATPipelines/PipelineRnaseq.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|ggplot2      |>=1.0.0     |1.0.0          |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|gplots       |>=2.14.2    |2.14.2         |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|grid         |>=3.1.1     |3.1.1          |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|gsnap        |>=2014-01-21|2014-01-21     |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|limma        |>=3.21.18   |3.21.18        |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|macs1        |>=1.4.2     |1.4.2          |False      |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|macs2        |>=2.0.10    |2.0.10.20131216|False      |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|peakranger   |>=1.16      |1.16           |False      |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|picardtools  |>=1.106     |1.106          |True       |CGATPipelines/PipelineMapping.py        |
|             |            |               |           |CGATPipelines/PipelineWindows.py        |
|             |            |               |           |CGATPipelines/pipeline_exome.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|reshape      |>=0.8.5     |0.8.5          |True       |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|sailfish     |>=0.6.3     |0.6.3          |True       |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|samr         |>=2.0       |2.0            |False      |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|samtools     |>=1.1       |1.1            |True       |CGATPipelines/PipelineMapping.py        |
|             |            |               |           |CGATPipelines/PipelinePeakcalling.py    |
|             |            |               |           |CGATPipelines/PipelineRnaseq.py         |
|             |            |               |           |CGATPipelines/PipelineWindows.py        |
|             |            |               |           |CGATPipelines/pipeline_exome.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|scripture    |>=2.0       |?              |True       |CGATPipelines/PipelinePeakcalling.py    |
+-------------+------------+---------------+-----------+----------------------------------------+
|siggenes     |>=1.39.0    |1.39.0         |False      |CGAT/Expression.py                      |
+-------------+------------+---------------+-----------+----------------------------------------+
|snow         |>=0.3.13    |0.3-13         |True       |scripts/runSPP.py                       |
+-------------+------------+---------------+-----------+----------------------------------------+
|snpEff       |>=4.0       |4.0e           |True       |CGATPipelines/pipeline_exome.py         |
+-------------+------------+---------------+-----------+----------------------------------------+
|spp          |>=?         |-              |True       |scripts/runSPP.py                       |
+-------------+------------+---------------+-----------+----------------------------------------+
|stampy       |>=1.0.23    |1.0.23         |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|star         |>=2.3.0e    |?              |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|tophat       |>=2.0.13    |2.0.13         |False      |CGATPipelines/PipelineMapping.py        |
+-------------+------------+---------------+-----------+----------------------------------------+
|ucsctools    |==?         |?              |True       |CGATPipelines/pipeline_windows.py       |
+-------------+------------+---------------+-----------+----------------------------------------+
|zinba        |>=2.01      |2.01           |True       |CGATPipelines/PipelinePeakcalling.py    |
|             |            |               |           |scripts/runZinba.py                     |
+-------------+------------+---------------+-----------+----------------------------------------+

What exactly is required will depend on the particular pipeline. The
pipeline assumes that the executables are in the users :envvar:`PATH`
and that the rest of the environment has been set up for each tool.

To check if the dependencies within a particular pipeline are satisfied, type::

   python cgat/CGATPipelines/pipeline_mapping.py check

To check all external dependencies, type::

   python scripts/cgat_list_dependencies.py

The dependencies are tracked through the module
:doc:`modules/Requirements`. Dependency tracking works by adding a
list of dependencies to the docstring of the module or script in which
the dependency arises.

Additionally, there is a list of additional software that is required
that are usually shipped as a source package with the operating system
such as sqlite_.

Please see the installation instructions for the `CGAT Toolkit
<http://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html>`_.

Python libraries
----------------

CGAT uses python extensively and is currently developed against python 2.7.1. Python
2.6 should work as well, but some libraries present in 2.7.1 but missing in 2.6
might need to be installed. Scripts have not yet been ported to python 3.

The CGAT pipelines require several python libraries to be installed.
When installing the CGAT code collection, these dependencies are
listed in the :file:` should be automatically installed

