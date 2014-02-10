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

Please see the installation instructions for the 
`CGAT Toolkit <http://www.cgat.org/~andreas/documentation/cgat/CGATInstallation.html>`_
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

The pipelines assume that Sun Grid Engine has been installed. Other queueing systems
might work, but expect to be disappointed. The pipeline is started on a 
:term:`submit host` assuming a default queue ``all.q``. Other queues can be specified on the
command line, for example::

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
versions. The following table gives a list software we have currently
(7/2/2014) installed:

.. How to create this table:
.. module -l list >& out; cat out | pe "s/ +/\t/g" | cut -f 1 | pe "s/\//\t/g" | tab2rst > x

+---------+---------------------+---------------+
|*Section*|*Software*           |*Version*      |
+---------+---------------------+---------------+
|apps     |java                 |jre1.6.0_26    |
+---------+---------------------+---------------+
|apps     |gccxml               |0.9            |
+---------+---------------------+---------------+
|apps     |R                    |2.15.2         |
+---------+---------------------+---------------+
|bio      |alignlib             |0.4.4          |
+---------+---------------------+---------------+
|apps     |python               |2.7.1          |
+---------+---------------------+---------------+
|apps     |perl                 |5.12.3         |
+---------+---------------------+---------------+
|apps     |graphlib             |0.1            |
+---------+---------------------+---------------+
|bio      |abacas               |1.3.1          |
+---------+---------------------+---------------+
|bio      |abiwtap              |1.2.1          |
+---------+---------------------+---------------+
|bio      |artemis              |15             |
+---------+---------------------+---------------+
|bio      |bamstats             |1.22           |
+---------+---------------------+---------------+
|bio      |bamtools1            |2.3.0          |
+---------+---------------------+---------------+
|bio      |batman               |0.2.3          |
+---------+---------------------+---------------+
|bio      |bedtools             |2.18.1         |
+---------+---------------------+---------------+
|bio      |belvu                |2.16           |
+---------+---------------------+---------------+
|bio      |bfast                |0.6.5a         |
+---------+---------------------+---------------+
|bio      |bioprospector        |2004           |
+---------+---------------------+---------------+
|bio      |boost                |1.46.1         |
+---------+---------------------+---------------+
|bio      |bowtie               |0.12.7         |
+---------+---------------------+---------------+
|bio      |bowtie2              |2.1.0          |
+---------+---------------------+---------------+
|bio      |BRIG-0.95-dist       |0.95           |
+---------+---------------------+---------------+
|bio      |BroadPeak            |1.0            |
+---------+---------------------+---------------+
|bio      |bwa                  |test-0.7.5a    |
+---------+---------------------+---------------+
|bio      |cdhit                |4.3            |
+---------+---------------------+---------------+
|bio      |CHANCE               |1.0            |
+---------+---------------------+---------------+
|bio      |circos               |0.63-4         |
+---------+---------------------+---------------+
|bio      |CLASS                |1.0.1          |
+---------+---------------------+---------------+
|bio      |clustalw             |2.1            |
+---------+---------------------+---------------+
|bio      |cortex               |0.0.4          |
+---------+---------------------+---------------+
|bio      |cortex_var           |1.0.5.20       |
+---------+---------------------+---------------+
|bio      |cufflinks            |2.0.2          |
+---------+---------------------+---------------+
|bio      |cpc                  |0.9-r2         |
+---------+---------------------+---------------+
|bio      |dialign              |2.2.1          |
+---------+---------------------+---------------+
|bio      |denovogear           |0.5.2          |
+---------+---------------------+---------------+
|bio      |emboss               |6.5.7          |
+---------+---------------------+---------------+
|bio      |ensembl              |62             |
+---------+---------------------+---------------+
|bio      |ensembl-variation    |62             |
+---------+---------------------+---------------+
|bio      |exonerate            |2.2.0          |
+---------+---------------------+---------------+
|bio      |famseq               |1.0.0-1        |
+---------+---------------------+---------------+
|bio      |fastqc               |0.9.2          |
+---------+---------------------+---------------+
|bio      |fastx                |0.0.13         |
+---------+---------------------+---------------+
|bio      |flash                |1.2.6          |
+---------+---------------------+---------------+
|bio      |flux-simulator       |1.2.1          |
+---------+---------------------+---------------+
|bio      |gatk                 |2.1.13         |
+---------+---------------------+---------------+
|bio      |gatk-full            |2.7-2          |
+---------+---------------------+---------------+
|bio      |gblocks              |0.91b          |
+---------+---------------------+---------------+
|bio      |gcprofile            |1.0            |
+---------+---------------------+---------------+
|bio      |gmap                 |2012.07.20     |
+---------+---------------------+---------------+
|bio      |galaxy               |dist           |
+---------+---------------------+---------------+
|bio      |gitools              |1.6.4          |
+---------+---------------------+---------------+
|bio      |glimmer              |3.02           |
+---------+---------------------+---------------+
|bio      |hpeak                |2.1            |
+---------+---------------------+---------------+
|bio      |gem                  |003425         |
+---------+---------------------+---------------+
|bio      |GEM                  |1.1            |
+---------+---------------------+---------------+
|bio      |GemSim               |version_unknown|
+---------+---------------------+---------------+
|bio      |idba                 |1.1.0          |
+---------+---------------------+---------------+
|bio      |IGV                  |2.3.2          |
+---------+---------------------+---------------+
|bio      |IGVTools             |2.1.24         |
+---------+---------------------+---------------+
|bio      |java_genomics_toolkit|0.0.1          |
+---------+---------------------+---------------+
|bio      |kent                 |1.0            |
+---------+---------------------+---------------+
|bio      |hmmer                |3.0            |
+---------+---------------------+---------------+
|bio      |lastz                |1.02.00        |
+---------+---------------------+---------------+
|bio      |leotools             |0.1            |
+---------+---------------------+---------------+
|bio      |mappability_map      |1.0            |
+---------+---------------------+---------------+
|bio      |Mauve                |2.3.1          |
+---------+---------------------+---------------+
|bio      |MEGAN                |4              |
+---------+---------------------+---------------+
|bio      |meme                 |4.9.1          |
+---------+---------------------+---------------+
|bio      |MetaGeneMark         |1.0.0          |
+---------+---------------------+---------------+
|bio      |metaphlan            |1.7.7          |
+---------+---------------------+---------------+
|bio      |meta-velvet          |1.2.02         |
+---------+---------------------+---------------+
|bio      |mtools               |1              |
+---------+---------------------+---------------+
|bio      |MUMmer               |3.23           |
+---------+---------------------+---------------+
|bio      |muscle               |3.8.31         |
+---------+---------------------+---------------+
|bio      |ncbiblast            |2.2.28+        |
+---------+---------------------+---------------+
|bio      |newickutils          |1.3.0          |
+---------+---------------------+---------------+
|bio      |novoalign            |2.07.11        |
+---------+---------------------+---------------+
|bio      |novoalignCS          |1.01.11        |
+---------+---------------------+---------------+
|bio      |paml                 |4.4c           |
+---------+---------------------+---------------+
|bio      |peakranger           |1.16           |
+---------+---------------------+---------------+
|bio      |peaksplitter         |1.0            |
+---------+---------------------+---------------+
|bio      |perm                 |0.3.5          |
+---------+---------------------+---------------+
|bio      |phylip               |3.69           |
+---------+---------------------+---------------+
|bio      |PhyloCSF             |1.0            |
+---------+---------------------+---------------+
|bio      |plinkseq             |0.08           |
+---------+---------------------+---------------+
|bio      |polymutt             |0.14           |
+---------+---------------------+---------------+
|bio      |polyphen             |2.0.23         |
+---------+---------------------+---------------+
|bio      |prodigal             |2.60           |
+---------+---------------------+---------------+
|bio      |Ray                  |2.2.0          |
+---------+---------------------+---------------+
|bio      |reaper               |13-100         |
+---------+---------------------+---------------+
|bio      |samtools             |0.1.19         |
+---------+---------------------+---------------+
|bio      |scripture            |2.0b           |
+---------+---------------------+---------------+
|bio      |seqimp               |13-095         |
+---------+---------------------+---------------+
|bio      |seqtk                |1.0.0          |
+---------+---------------------+---------------+
|bio      |sga                  |1.0.0          |
+---------+---------------------+---------------+
|bio      |shrimp               |2.1.1          |
+---------+---------------------+---------------+
|bio      |sicer                |1.1            |
+---------+---------------------+---------------+
|bio      |sickle               |1.0            |
+---------+---------------------+---------------+
|bio      |sift                 |4.0.3          |
+---------+---------------------+---------------+
|bio      |simseq               |72ce499        |
+---------+---------------------+---------------+
|bio      |snpEff               |3.3            |
+---------+---------------------+---------------+
|bio      |soap                 |2.21           |
+---------+---------------------+---------------+
|bio      |SOAPdenovo2          |2.04           |
+---------+---------------------+---------------+
|bio      |soapsplice           |1.0            |
+---------+---------------------+---------------+
|bio      |SPAdes-3.0           |3.0            |
+---------+---------------------+---------------+
|bio      |sratoolkit           |2.1.7          |
+---------+---------------------+---------------+
|bio      |SpliceMap            |3.3.5.2        |
+---------+---------------------+---------------+
|bio      |stampy               |1.0.17         |
+---------+---------------------+---------------+
|bio      |statgen              |0.1.4          |
+---------+---------------------+---------------+
|bio      |star                 |2.3.0e         |
+---------+---------------------+---------------+
|bio      |storm                |0.1            |
+---------+---------------------+---------------+
|bio      |subread              |1.3.6-p1       |
+---------+---------------------+---------------+
|bio      |sylamer              |08-123         |
+---------+---------------------+---------------+
|bio      |tabix                |0.2.6          |
+---------+---------------------+---------------+
|bio      |TAXAassign           |0.4            |
+---------+---------------------+---------------+
|bio      |triodenovo           |0.01           |
+---------+---------------------+---------------+
|bio      |tophat               |1.4.1          |
+---------+---------------------+---------------+
|bio      |tophat2              |2.0.10         |
+---------+---------------------+---------------+
|bio      |treebest             |0.1            |
+---------+---------------------+---------------+
|bio      |Trinity              |2012-01-25     |
+---------+---------------------+---------------+
|bio      |tv                   |0.5            |
+---------+---------------------+---------------+
|bio      |vcftools             |0.1.8a         |
+---------+---------------------+---------------+
|bio      |velvet               |1.2.10         |
+---------+---------------------+---------------+
|bio      |velvet-optimiser     |2.2.5          |
+---------+---------------------+---------------+
|bio      |VEP                  |67             |
+---------+---------------------+---------------+

What exactly is required will depend on the particular pipeline. The pipeline assumes
that the executables are in the users :envvar:`PATH` and that the rest of the environment
has been set up for each tool.

Additionally, there is a list of additional software that is required
that are usually shipped as a source package with the operating
system such as sqlite_.

Please see the installation instructions for the 
`CGAT Toolkit <http://www.cgat.org/~andreas/documentation/cgat/CGATInstallation.html>`_.

Python libraries
----------------

CGAT uses python extensively and is currently developed against python 2.7.1. Python
2.6 should work as well, but some libraries present in 2.7.1 but missing in 2.6
might need to be installed. Scripts have not yet been ported to python 3.

The CGAT pipelines require several python libraries to be installed.
When installing the CGAT code collection, these dependencies should be
automatically installed.

