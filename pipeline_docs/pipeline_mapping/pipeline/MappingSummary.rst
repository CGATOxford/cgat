===============
Mapping summary
===============

This section gives an overview over the mapping process. 

Alignment summary
=================

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.BamSummary
   :render: table
   :slices: mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :slices: mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per number of mismatches in alignment.

Reads summary
=============

The following table 

.. report:: Mapping.BamSummary
   :render: table
   :slices: reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :slices: reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.MappingFlagsHits
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of reads per number of alignments (hits) per read.

