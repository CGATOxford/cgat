*****************************************
GPipe - Gene prediction pipeline
*****************************************

Introduction
------------

This document describes the pipeline of the Chris Ponting group
for predicting genes by homology. The input is a set of known transcripts
from a reference genome and the masked genomic sequence of a
target genome.

The pipeline predicts genes in a two-step procedure:

1. Regions of similarity between the known transcripts and the reference genome 
   are identified using a quick heuristic search.

2. The regions of similarity of step 1 are submitted to a sensitive, but slower gene 
   prediction program.

The pipeline contains many options to mask sequences, analyse and quality
control the predictions and store the results in a relational database.
This manual describes how to setup up the pipeline, run it on our cluster,
and analyse the results.

.. note::
   The pipeline is tailored to the computational setup within the Ponting group.
   Porting it elsewhere will require a significant amount of software installation
   and configuration. 

Setting up
----------

The pipeline consists of a set of scripts and a makefile,
that glues together the various scripts. This section tells
you what programs are needed to be installed (`Requirements`_),
how to install the software from the pipeline (`Installation`_).

Next, the input files need to prepared (`Input files`_).

Before starting the pipeline, you need to configure your environment
and the pipeline (`Configuration`_).

Requirements
++++++++++++

Gpipe requires following software to be installed:

   postgres
       Gpipe currently requires a postgres installation.

   seq
       Seg is a program to mask low complexity regions in protein sequences.

   exonerate
       Exonerate is a program to align a peptide sequence to a genomic sequence 
       (other alignment modes are possible). It offers heuristic modes, that allow for fast 
       scanning of large chunks of genomic DNA, and exhaustive modes, that do a full dynamic
       programming mode. It is available at http://www.ebi.ac.uk/~guy/exonerate.

   python
       Most scripts require python...

   perl
      ...unless they are written in perl

   alignlib
      a library for sequence alignments and its python interface. See
      http://sourceforge.net/projects/alignlib.

   SGE
      Sun Grid Engine. Other schedulers might work.

   seq_pairs_kaks
      Leo's wrapper around PAML (optional, only required for step12).

Installation
++++++++++++

Untar and unpack the source code tarball. The directory in which it ends up
is the :term:`source directory`.


Gpipe works in a :term:`working directory`. To set up the pipeline with the current 
directory as the :term:`working directory`, run:: 

   python <src>setup.py --method=gpipe --project=project_name > setup.log

:file:`src` is the location of the :term:`source directory`.

This will create a makefile in the current directory and give the project the
name ``project_name``. The latter will be used as the name of the database schema
in postgres in which data will be stored.

Note that gpipe assumes that you use ``bash`` as you shell. In particular, it requires
that two functions are in your environment::

   #-----------------------------------
   # helper functions for detecting errors in pipes
   #-----------------------------------
   detect_pipe_error_helper()
   {
    while [ "$#" != 0 ] ; do
        # there was an error in at least one program of the pipe
        if [ "$1" != 0 ] ; then return 1 ; fi
        shift 1
    done
    return 0
   }

   detect_pipe_error()
   { 
    detect_pipe_error_helper "${PIPESTATUS[@]}"
    return $?
   }

Once the code is in place, add the input files to the working and make sure
that all the other requirements are fulfilled.

Input files
+++++++++++

Gpipe requires 5 files to run. These input files contain the reference gene set to predict with and
the genome sequence to predict in. Sample data is available at 
http://genserv.anat.ox.ac.uk/downloads/software/gpipe/sampledata/gpipe_sample_data.tar.
To drop the sample data into your :term:`working directory`, type::

   wget ttp://genserv.anat.ox.ac.uk/downloads/software/gpipe/sampledata/gpipe_sample_data.tar
   tar -xf gpipe_sample_data.tar
   gunzip *

The filenames and their contents are:

:file:`peptides.fasta`
   A fasta-formatted file with peptide sequences. Each sequence is on a single line. 
   The identifier of a sequence is taken from the description line 
   with the pattern ``>(\S+)`` (characters between > and first white-space).
   For example::
   
      >CG11023-RA
      MGERDQPQSSERISIFNPPVYTQHQVRNEAPYIPTTFDLLSDDEESSQRVANAGPSFRPL...
      >CG2671-RA
      MLKFIRGKGQQPSADRHRLQKDLFAYRKTAQHGFPHKPSALAYDPVLKLMAIGTQTGALK...
      ...

:file:`genome.fasta` and :file:`genome.idx`
   A fasta-formatted file with the genomic sequence together with its index. This file can be created
   from a collection of individual fasta formatted files (for example, per chromosome files) using the command::
   
      python <src>index_fasta.py genome <somedir>my_dna_sequences*.fa.gz > genome.log
      
:file:`reference.exons`
   A table with gene models from the reference gene set. 
   This is a tab-formatted table with the following columns:

   transcript name	
      name of the transcript consistent with :file:`peptides.fasta`
   contig name
      name of the DNA segment the transcript is located on
   strand
      strand
   phase
      phase of that particular exon
   exon-id
      numerical number of exon
   peptide_start
      start of exon in transcript sequence 
   peptide_end
      end of exon in transcript sequence
   genome_start
      start of exon on contig
   genome_end
      end of exon on contig

   Coordinates are 0-based, half-open intervals. Genomic coordinates are forward/reverse strand coordinates.

   For example::

     CG10000-RA	chr3R	-1	0	1	0	126	24577165	24577291
     CG10000-RA	chr3R	-1	0	2	126	287	24576946	24577107
     CG10000-RA	chr3R	-1	1	3	287	466	24576706	24576885
     CG10000-RA	chr3R	-1	2	4	466	930	24576187	24576651
     CG10000-RA	chr3R	-1	0	5	930	1100	24575892	24576062
     CG10000-RA	chr3R	-1	1	6	1100	1280	24575573	24575753
     CG10000-RA	chr3R	-1	1	7	1280	1677	24574936	24575330
     CG10001-RA	chr3R	-1	0	1	0	540	24569207	24569747
     CG10001-RA	chr3R	-1	0	2	540	819	24566427	24566706
     CG10001-RA	chr3R	-1	0	3	819	978	24566193	24566352

:file:`map_rep2mem`
   A table linking genes to transcripts. This tab-formatted table contains the following columns
  
   rep
      A gene identifier
   mem
      A transcript identifier
   size
      Transcript size

Configuration
+++++++++++++

To configure the pipeline, options can be set in the :file:`Makefile` in the
:term:`working directory`.

Options that might need to be changed:

   PARAM_PSQL_DATABASE
      The psql database

   PARAM_PSQL_HOST
      The psql host

   PARAM_PSQL_USER
      The psql username

Running the pipeline
--------------------

The pipeline uses makefiles to control script logic. Before executing any make commands,
run::

        source setup.csh

to update your paths and other environment variables.

Before first running the pipeline, some maintenance work need to
be performed like creating the database schema and the tables.
To prepare the pipeline, run::

   make prepare

This needs to be done only once. 

To run the pipeline, type::

   make all

Gpipe writes status messages to the file :file:`log`
in the :term:`working directory`.

Results
-------

Results of the gpipe run are stored in the psql database.

The view ``overview`` aggregates most results into a single
table for easy access.

Troubleshooting
---------------

If something goes wrong, the first step is to look at
the command line that caused the problem. To see the command
executed, run::

   make -n <target>

We use Sun Grid Engine as job queueing system and assume that for all nodes
the code and data can be reached via the same mount point. All jobs that are run on
the cluster are prefixed by the MAKE variable $(CMD_REMOTE_SUBMIT). You can set this
variable to the empty string to run everything locally or on a mosix cluster:

        make all CMD_REMOTE_SUBMIT=


Steps
-----

The pipeline proceeds in 12 steps, which are:

Step1
	Masking of protein sequences  
Step2
	Selecting representative transcipts to search with
Step3
	Running exonerate    
Step4
	Running TBLASTN (disabled)
Step5	
	Collating putative gene-containing regions
Step6	
	Predicting genes for representative transcripts.
Step7
	Predicting genes for redundant (alternative) transcripts
Step8	
	Predicting genes for member sequences  
Step9	
	Analysing the predictions  
Step10
	Quality control of predictions  
Step11
        Removing redundant/erroneous predictions
Step12
        Filter by ks (optional)

Glossary
---------
.. glossary::

   working directory
     The working directory. Location of the data files and results. All commands in this tutorial are 
     executed in the working directory.

   source directory
     The location of the source code. The place where the script :file:`setup.py` resides.
