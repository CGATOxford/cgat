.. _CGATSetup:

=========================
Installing CGAT pipelines
=========================

The CGAT pipelines, scripts and libraries make several assumptions about
the computing environment. This section describes how to install the code
and set up your computing environment.

Downloading and installing the source code
==========================================

To obtain the latest code, check it out from the public mercurial_ repository at::

   hg clone http://www.cgat.org/hg/cgat/ src

Once checked-out, you can get the latest changes via pulling and updating::

   hg pull 
   hg update

Some scripts contain cython code that needs to be recompiled if the
script or the pysam_ installation has changed. To rebuild all scripts,
for example after updating the repository, type::

   python src/rebuild_extensions.py

Recompilation requires a C compiler to be installed. 

Setting up the computing environment
====================================

The pipelines assume that Sun Grid Engine has been installed. Other queueing systems
might work, but expect to be disappointed. The pipeline is started on a 
:term:`submit host` assuming a default queue ``all.q``. Other queues can be specified on the
command line, for example::

    python <src>pipeline_<name>.py --cluster-queue=medium_jobs.q

A pipeline might start up to ``-p/--multiprocess`` processes. Preferentially,
tasks are sent to the cluster, but for some tasks this is not possible. 
These might thus run on the :term:`submit host`, so make sure it is fairly powerful.

Pipelines expects that the :term:`working directory` is accessible with
the same path both from the submit and the :term:`execution host`. 

Software requirements
=====================

On top of pipeline specific bioinformatics software, CGAT pipelines make use of the 
following software:

+----------------------+-------------------+------------------------------------------+
|*Program*             |*Version*          |*Purpose*                                 |
+----------------------+-------------------+------------------------------------------+
|bowtie_               |>=0.12.7           |read mapping                              |
+----------------------+-------------------+------------------------------------------+
|tophat_               |>=1.3.0            |read mapping                              |
+----------------------+-------------------+------------------------------------------+
|cufflinks_            |>=0.9.3            |transcription levels                      |
+----------------------+-------------------+------------------------------------------+
|samtools              |>=0.1.12           |bam/sam files                             |
+----------------------+-------------------+------------------------------------------+
|Picard                |                   |bam/sam files                             |
+----------------------+-------------------+------------------------------------------+
|bedtools              |                   |working with intervals                    |
+----------------------+-------------------+------------------------------------------+
|sra-tools             |                   |extracting reads from .sra files          |
+----------------------+-------------------+------------------------------------------+
|Jim Kent's UCSC tools |                   |working with genomic data                 |
+----------------------+-------------------+------------------------------------------+
|R + bioconductor      |                   |statistical functions.                    |
+----------------------+-------------------+------------------------------------------+
|                      |                   |                                          |
+----------------------+-------------------+------------------------------------------+

What exactly is required will depend on the particular pipeline. The pipeline assumes
that the executables are in the users :envvar:`PATH` and that the rest of the environment
has been set up for each tool.

Python libraries
----------------

CGAT uses python extensively and is currently developed agains python 2.7.1. Python
2.6 should work as well, but some libraries present in 2.7.1 but missing in 2.6
might need to be installed. Scripts have not yet been ported to python 3.

CGAT requires the following in-house python libraries to be installed:

+--------------------+-------------------+----------------------------------------+
|*Library*           |*Version*          |*Purpose*                               |
+--------------------+-------------------+----------------------------------------+
|pysam_              |0.5.0              |python bindings for samtools            |
+--------------------+-------------------+----------------------------------------+
|NCL_                |0.1                |nested containment lists                |
+--------------------+-------------------+----------------------------------------+
|fastgtf_            |0.1                |fast gtf parsing                        |
+--------------------+-------------------+----------------------------------------+
|alignlib_           |0.4.4              |C++ sequence alignment library with     |
|                    |                   |python bindings.                        |
+--------------------+-------------------+----------------------------------------+
|Components          |0.1                |connected components computation        |
+--------------------+-------------------+----------------------------------------+
|sphinxreport_       |latest             |report generator                        |
+--------------------+-------------------+----------------------------------------+

In addition, CGAT scripts make extensive use of the following python libraries (list below
might not be complete):

+--------------------+-------------------+----------------------------------------+
|*Library*           |*Version*          |*Purpose*                               |
+--------------------+-------------------+----------------------------------------+
|numpy               |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|scipy               |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|rpy2                |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|matplotlib          |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|ruffus              |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|drmaa_python        |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|                    |                   |                                        |
+--------------------+-------------------+----------------------------------------+
|                    |                   |                                        |
+--------------------+-------------------+----------------------------------------+

.. _alignlib: http://wwwfgu.anat.ox.ac.uk/~andreas/alignlib
.. _ncl: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/ncl/contents.html
.. _fastgtf: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/fastgtf/contents.html
.. _pysam: http://code.google.com/p/pysam/
.. _sphinxreport: http://code.google.com/p/sphinx-report/
.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _mercurial: http://mercurial.selenic.com
