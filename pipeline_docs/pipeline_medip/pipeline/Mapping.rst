============
Read Mapping
============

The following table presents an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,reverse,duplicate

   Mapping Summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reverse,duplicate

   Mapping Summary

.. report:: Mapping.PicardAlign
   :render: table
   :force:

   Picard Alignment Statistics

.. report:: Mapping.PicardAlignPlot
   :render: interleaved-bar-plot

   Picard Alignment Statistics

