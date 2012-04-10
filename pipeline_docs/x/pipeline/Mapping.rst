=======
Mapping
=======

Mapping results
===============


The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,unmapped,reverse,duplicate

   Mapping Summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,unmapped,reverse,duplicate

   Mapping Summary

.. report:: Mapping.PicardAlign
   :render: table
   :tracks: picard_align_stats
   :slices: Total_read_pairs,Aligned_pairs,Strand_balance,Duplicates

   Picard Alignment Statistics

.. report:: Mapping.PicardAlignPlot
   :render: interleaved-bar-plot
   :slices: Total_read_pairs
   :width:  1000

   Picard Alignment Statistics

