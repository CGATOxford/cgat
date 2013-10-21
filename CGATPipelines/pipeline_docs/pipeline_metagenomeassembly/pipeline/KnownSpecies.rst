.. _knownSpecies:


======================================
Mapping of reads to reference genome
======================================


From taxonomic profiling of raw reads, we get a rough idea of the species that are present. Many of these
are not classified at the level of the species which is consistent with diversification of the sample
from reference genetic markers. It is not easy to assemble complete or near-complete genomes from 
metagenomic data and so as a first-pass it is often desirable to assess how many reads align to genomes
that are present in a reference genome database such as the NCBI. To ensure that we do not waste too 
much time on mapping to the thousands of reference genomes that are currently available, we utilise
the information derived from taxonomic profiles. All reference genomes - and if desirable plasmids - that
are from any of the genera that are identified in the taxonomic profiling analysis are used for alignment.
Here we describe the genomes that are represented and how many alignments we see.



Reference genomes
==================

The reference genomes used are described below.


.. report:: knownSpecies.SpeciesCount
   :render: table
   

   Total number of reference genomes aligned to


.. report:: knownSpecies.Species
   :render: table
   

   List of reference genomes aligned to




Alignment summary to reference genomes
=======================================

Below is a summary of the results of alignments to reference genomes


.. report:: knownSpecies.KnownAlignments
   :render: interleaved-bar-plot
   

   Alignment summary for reference genome alignment






