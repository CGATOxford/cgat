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

On top of pipeline specific bioinformatics software, CGAT pipelines
make use a variety of software. Unfortunately we can't support many
versions. The following table gives a list software we have currently
installed:

.. How to create this table:
.. module -l list >& out; cat out | pe "s/ +/\t/g" | cut -f 1 | pe "s/\//\t/g" | tab2rst > x

+---------+-----------------+-----------+
|*Section*|*Software*       |*Version*  |
+---------+-----------------+-----------+
|apps     |java             |jre1.6.0_26|
+---------+-----------------+-----------+
|apps     |gccxml           |0.9        |
+---------+-----------------+-----------+
|apps     |R                |2.14.1     |
+---------+-----------------+-----------+
|bio      |alignlib         |0.4.4      |
+---------+-----------------+-----------+
|apps     |python           |2.7.1      |
+---------+-----------------+-----------+
|apps     |perl             |5.12.3     |
+---------+-----------------+-----------+
|apps     |graphlib         |0.1        |
+---------+-----------------+-----------+
|bio      |abiwtap          |1.2.1      |
+---------+-----------------+-----------+
|bio      |bamstats         |1.22       |
+---------+-----------------+-----------+
|bio      |batman           |0.2.3      |
+---------+-----------------+-----------+
|bio      |bedtools         |2.13.3     |
+---------+-----------------+-----------+
|bio      |belvu            |2.16       |
+---------+-----------------+-----------+
|bio      |bfast            |0.6.5a     |
+---------+-----------------+-----------+
|bio      |bioprospector    |2004       |
+---------+-----------------+-----------+
|bio      |bowtie           |0.12.7     |
+---------+-----------------+-----------+
|bio      |bwa              |0.5.9      |
+---------+-----------------+-----------+
|bio      |cdhit            |4.3        |
+---------+-----------------+-----------+
|bio      |clustalw         |2.1        |
+---------+-----------------+-----------+
|bio      |cufflinks        |1.3.0      |
+---------+-----------------+-----------+
|bio      |cpc              |0.9-r2     |
+---------+-----------------+-----------+
|bio      |dialign          |2.2.1      |
+---------+-----------------+-----------+
|bio      |ensembl          |62         |
+---------+-----------------+-----------+
|bio      |ensembl-variation|62         |
+---------+-----------------+-----------+
|bio      |exonerate        |2.2.0      |
+---------+-----------------+-----------+
|bio      |fastqc           |0.9.2      |
+---------+-----------------+-----------+
|bio      |fastx            |0.0.13     |
+---------+-----------------+-----------+
|bio      |gatk             |1.0.5506   |
+---------+-----------------+-----------+
|bio      |gblocks          |0.91b      |
+---------+-----------------+-----------+
|bio      |gcprofile        |1.0        |
+---------+-----------------+-----------+
|bio      |gmap             |2011.03.28 |
+---------+-----------------+-----------+
|bio      |galaxy           |dist       |
+---------+-----------------+-----------+
|bio      |IGV              |2.0.23     |
+---------+-----------------+-----------+
|bio      |IGVTools         |1.5.12     |
+---------+-----------------+-----------+
|bio      |kent             |1.0        |
+---------+-----------------+-----------+
|bio      |hmmer            |3.0        |
+---------+-----------------+-----------+
|bio      |leotools         |0.1        |
+---------+-----------------+-----------+
|bio      |meme             |4.7.0      |
+---------+-----------------+-----------+
|bio      |muscle           |3.8.31     |
+---------+-----------------+-----------+
|bio      |mappability_map  |1.0        |
+---------+-----------------+-----------+
|bio      |ncbiblast        |2.2.25+    |
+---------+-----------------+-----------+
|bio      |newickutils      |1.3.0      |
+---------+-----------------+-----------+
|bio      |novoalign        |2.07.11    |
+---------+-----------------+-----------+
|bio      |novoalignCS      |1.01.11    |
+---------+-----------------+-----------+
|bio      |paml             |4.4c       |
+---------+-----------------+-----------+
|bio      |picard-tools     |1.48       |
+---------+-----------------+-----------+
|bio      |phylip           |3.69       |
+---------+-----------------+-----------+
|bio      |polyphen         |2.0.23     |
+---------+-----------------+-----------+
|bio      |samtools         |0.1.18     |
+---------+-----------------+-----------+
|bio      |shrimp           |2.1.1      |
+---------+-----------------+-----------+
|bio      |sicer            |1.1        |
+---------+-----------------+-----------+
|bio      |sift             |4.0.3      |
+---------+-----------------+-----------+
|bio      |simseq           |72ce499    |
+---------+-----------------+-----------+
|bio      |soap             |2.21       |
+---------+-----------------+-----------+
|bio      |soapsplice       |1.0        |
+---------+-----------------+-----------+
|bio      |sratoolkit       |2.1.7      |
+---------+-----------------+-----------+
|bio      |SpliceMap        |3.3.5.2    |
+---------+-----------------+-----------+
|bio      |stampy           |1.0.17     |
+---------+-----------------+-----------+
|bio      |statgen          |0.1.4      |
+---------+-----------------+-----------+
|bio      |storm            |0.1        |
+---------+-----------------+-----------+
|bio      |tabix            |0.2.5      |
+---------+-----------------+-----------+
|bio      |tophat           |1.4.1      |
+---------+-----------------+-----------+
|bio      |treebest         |0.1        |
+---------+-----------------+-----------+
|bio      |tv               |0.5        |
+---------+-----------------+-----------+
|bio      |vcftools         |0.1.8a     |
+---------+-----------------+-----------+
|bio      |emboss           |6.3.1      |
+---------+-----------------+-----------+
|bio      |velvet           |1.1.04     |
+---------+-----------------+-----------+
|bio      |perm             |0.3.5      |
+---------+-----------------+-----------+
|bio      |lastz            |1.02.00    |
+---------+-----------------+-----------+
|bio      |hpeak            |2.1        |
+---------+-----------------+-----------+
|bio      |boost            |1.46.1     |
+---------+-----------------+-----------+
|bio      |Trinity          |2012-01-25 |
+---------+-----------------+-----------+
|bio      |bowtie2          |2.0.0-beta5|
+---------+-----------------+-----------+
|bio      |tophat2          |2.0.0      |
+---------+-----------------+-----------+
|bio      |all              |1.0        |
+---------+-----------------+-----------+

What exactly is required will depend on the particular pipeline. The pipeline assumes
that the executables are in the users :envvar:`PATH` and that the rest of the environment
has been set up for each tool.

Additionally, there is a list of additional software that is required
that are usually shipped as a source package with the operating
system. These are:

sqlite

Python libraries
----------------

CGAT uses python extensively and is currently developed against python 2.7.1. Python
2.6 should work as well, but some libraries present in 2.7.1 but missing in 2.6
might need to be installed. Scripts have not yet been ported to python 3.

CGAT requires the following in-house python libraries to be installed:

+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|*Library*           |*Version*          |*Purpose*                               |*Download*                                                                                                                     |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|pysam_              |0.6.0              |python bindings for samtools            |hg clone https://code.google.com/p/pysam/ pysam                                                                                |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|ncl_                |0.1                |nested containment lists                |hg clone http://www.cgat.org/hg/ncl ncl                                                                                        |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|fastgtf_            |0.1                |fast gtf parsing                        |hg clone http:://www.cgat.org/hg/fastgtf fastgtf                                                                               |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|alignlib_           |0.4.5              |C++ sequence alignment library with     |wget http://downloads.sourceforge.net/project/alignlib/alignlib/alignlib-0.4.5.tar.gz                                          |
|                    |                   |python bindings.                        |                                                                                                                               |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|components          |0.1                |connected components computation        |hg clone http://www.cgat.org/hg/components components                                                                          |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+
|sphinxreport_       |latest             |report generator                        |svn checkout https://sphinx-report.googlecode.com/svn/trunk/ sphinx-report                                                     |
+--------------------+-------------------+----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+

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

The full list of modules installed at CGAT is:

+-------------+--------+-------+
|Module       |Version |Method |
+-------------+--------+-------+
|pycairo      |01/08/06|S      |
+-------------+--------+-------+
|pygjobject   |2.20.0  |S      |
+-------------+--------+-------+
|pygtk        |2.16.0  |S      |
+-------------+--------+-------+
|wxPython     |2.9.1.1 |S      |
+-------------+--------+-------+
|matplotlib   |1       |S      |
+-------------+--------+-------+
|numpy        |01/05/01|E      |
+-------------+--------+-------+
|scipy        |0.8.0   |S      |
+-------------+--------+-------+
|rpy          |1.0.3   |S      |
+-------------+--------+-------+
|rpy2         |02/02/00|S      |
+-------------+--------+-------+
|networkx     |1.3     |E      |
+-------------+--------+-------+
|pytables     |2.2     |       |
+-------------+--------+-------+
|pygccxml     |1       |S      |
+-------------+--------+-------+
|pyplusplus   |1       |S      |
+-------------+--------+-------+
|bx.python    |        |       |
+-------------+--------+-------+
|pygresql     |4       |E      |
+-------------+--------+-------+
|myqsl-python |01/02/03|E      |
+-------------+--------+-------+
|biopython    |1.56    |E      |
+-------------+--------+-------+
|ply          |3.3     |E      |
+-------------+--------+-------+
|psyco        |        |       |
+-------------+--------+-------+
|pyrex        |0.9.9   |E      |
+-------------+--------+-------+
|cython       |0.13    |E      |
+-------------+--------+-------+
|sphinx       |1.0.5   |E      |
+-------------+--------+-------+
|reportlab    |2.5     |E      |
+-------------+--------+-------+
|guppy        |0.1.9   |E      |
+-------------+--------+-------+
|pil          |01/01/07|E      |
+-------------+--------+-------+
|threadpool   |01/02/07|E      |
+-------------+--------+-------+
|progressbar  |2.3     |E      |
+-------------+--------+-------+
|virtualenv   |01/05/01|E      |
+-------------+--------+-------+
|sqlalchemy   |0.6.5   |E      |
+-------------+--------+-------+
|ruffus       |2.2     |E      |
+-------------+--------+-------+
|drmaa        |0.4b3   |E      |
+-------------+--------+-------+
|bx.python    |12/01/10|S      |
+-------------+--------+-------+
|corebio      |0.5.0   |E      |
+-------------+--------+-------+
|weblogolib   |3       |E      |
+-------------+--------+-------+
|NCL          |        |C      |
+-------------+--------+-------+
|fastgtf      |        |C      |
+-------------+--------+-------+
|Components   |        |C      |
+-------------+--------+-------+
|mercurial    |01/07/03|E      |
+-------------+--------+-------+
|scikits.learn|0.7.1   |E      |
+-------------+--------+-------+
|web.py       |0.34    |E      |
+-------------+--------+-------+
|pandas       |0.5.0   |E      |
+-------------+--------+-------+
|pybedtools   |0.6     |E      |
+-------------+--------+-------+

Method : Installation method (E = easy_install/setuptools, S =
setup.py/distutils, C = CGAT)

.. _alignlib: http://wwwfgu.anat.ox.ac.uk/~andreas/alignlib
.. _ncl: http://www.cgat.org/~andreas/documentation/ncl/contents.html
.. _fastgtf: http://www.cgat.org/~andreas/documentation/fastgtf/contents.html
.. _pysam: http://code.google.com/p/pysam/
.. _sphinxreport: http://code.google.com/p/sphinx-report/
.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _mercurial: http://mercurial.selenic.com
