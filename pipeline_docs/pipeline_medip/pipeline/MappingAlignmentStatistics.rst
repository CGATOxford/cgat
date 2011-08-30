====================
Alignment statistics
====================

This section brings together various alignment statistics.

Alignments
==========

The following table present an overview of the alignments in the 
BAM files for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.PicardAlign
   :render: table
   :force:

   Picard Alignment Statistics

.. report:: Mapping.PicardAlignPlot
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 200

   Picard Alignment Statistics

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :tracks: r(.genome)
   :render: table
   :force:

   Alignments summary

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :tracks: r(.genome)
   :render: interleaved-bar-plot
   :slices: PCT_PF_READS_ALIGNED,STRAND_BALANCE

   Percentage quantities

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :tracks: r(.genome)
   :render: interleaved-bar-plot
   :slices: PF_READS_ALIGNED,PF_HQ_ALIGNED_READS

   Percentage quantities

Quality
=======

.. report:: Mapping.PicardQualityByCycleHistogram
   :tracks: r(.genome)
   :render: line-plot
   :as-lines:
   :yrange: 0,

   mean quality score by cycle

.. report:: Mapping.PicardQualityDistributionHistogram
   :tracks: r(.genome)
   :render: line-plot
   :as-lines:
   :yrange: 0,

   quality score distribution

Duplicates
===========

.. report:: Mapping.PicardDuplicatesMetrics
   :tracks: r(.prep)
   :render: table
   :force:

   Insert size summary

.. report:: Mapping.PicardDuplicatesHistogram
   :tracks: r(.prep)
   :render: line-plot
   :as-lines:
   :yrange: 0,

   Histogram of insert sizes

Insert sizes
============

.. report:: Mapping.PicardInsertSizeMetrics
   :render: table
   :force:

   Insert size summary

.. report:: Mapping.PicardInsertSizeHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,

   Histogram of insert sizes
