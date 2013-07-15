*******************************
454 Transcript mapping pipeline
*******************************

Purpose
-------

Map 454 reads onto a genome and assemble overlapping
transcripts into transcript models.

The pipeline currently does not use base quality information
during mapping and does not consider alternative transcripts.


Setting up
----------

To set up the pipeline in the current directory run::

   python setup.py --method=map_transcripts_454 > setup.log

Add or link fasta files of reads into directory. These should end
with the suffix .fasta. The pipeline will process several files at the same time.
For example::

   tissue1.fasta
   tissue2.fasta
   tissue3.fasta

Link towards the genome from /net/cpp-data/backup/databases/indexed_fasta and
call the files genome.fasta and genome.idx. For example::
       
   ln -s /net/cpp-data/backup/databases/indexed_fasta/hs_ncbi36_softmasked.fasta genome.fasta
   ln -s /net/cpp-data/backup/databases/indexed_fasta/hs_ncbi36_softmasked.idx genome.idx

Build the index for `gmap`_ by running gmap_setup. By default, `gmap`_ indices should be put
in :file:`/net/cpp-mirror/databases/gmap`. Provide the location to the indices using
the variable ``PARAM_GMAP_OPTIONS``.

.. note::

   Indices on networked disks are slow to load up. For performance reasons 
   work with local indices.

Configuration
-------------

Edit the :file:`Makefile` to configure the pipeline. See Parameters_ below.

Parameters
----------

The following parameters can be set in the :file:`Makefile`:

.. report:: Trackers.MakefileParameters
   :render: glossary
   :tracks: Makefile.map_transcripts_454
   :transpose:

   Overview of pipeline parameters.


.. _gmap: http://www.molecularevolution.org/software/genomics/gmap
