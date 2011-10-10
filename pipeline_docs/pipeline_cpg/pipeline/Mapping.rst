===================
Bowtie Read Mapping
===================

Bam2stats
---------

The following table presents an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,reverse,duplicates,unique_alignments

   Bam2stats Mapping Statistics

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reverse,duplicates,unique_alignments

   Bam2stats Mapping Statistics

Picard
------

.. report:: Mapping.PicardAlignPerc
   :render: table
   :force:

   Picard Alignment Percentages

.. report:: Mapping.PicardAlignCount
   :render: table
   :force:

   Picard Alignment Statistics

.. report:: Mapping.PicardAlignPlot
   :render: interleaved-bar-plot

   Picard Alignment Statistics

