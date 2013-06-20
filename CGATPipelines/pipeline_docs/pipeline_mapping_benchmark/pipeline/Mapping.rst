=======
Mapping
=======

Mapping results
===============

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,unmapped,reads_mapped,reverse,reads_unique_alignments,duplicate,total,mapped

   Mapping Summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reverse,reads_unique_alignments,duplicate
   :width:  1200
   :force:

   Mapping Summary

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per number of mismatches in alignment.

.. report:: Mapping.MappingFlagsHits
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per read.


