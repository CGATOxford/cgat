==================
Mapping statistics
==================

Alignment metrics
==============

The following table present an overview of the alignments in the 
BAM files for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :render: xls-table

   Alignments summary

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :render: interleaved-bar-plot
   :slices: PCT_PF_READS_ALIGNED,STRAND_BALANCE
   :layout: row
   :split-at: 10

   Percentage quantities

.. report:: Mapping.PicardAlignmentSummaryMetrics
   :render: interleaved-bar-plot
   :slices: PF_READS_ALIGNED,PF_HQ_ALIGNED_READS
   :layout: row
   :split-at: 10

   Percentage quantities

Quality metrics
===============

The following plots show the distribution of base quality scores in
reads. Note that if the bam-files have been stripped of sequence and
quality information the plots below will contain meaningless values.

.. report:: Mapping.PicardQualityByCycleHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :layout: row

   mean quality score by cycle

.. report:: Mapping.PicardQualityDistributionHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :layout: row

   quality score distribution

Insert size metrics
===================

The following table presents an overview of the insert size metrics
for each :term:`track`.  See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#InsertSizeMetrics>`_
for a definition of the field contents.
(note: only applicable for paired-end data).

.. report:: Mapping.PicardInsertSizeMetrics
   :render: xls-table

   Insert size summary

.. report:: Mapping.PicardInsertSizeHistogram
   :render: line-plot
   :as-lines:

   Insert size distribution

Duplicate metrics
=================

The current section is currently not completed. No data will be shown below.


The following table presents an overview of the duplicate metrics
for each :term:`track`.  See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#DuplicationMetrics>`_
for a definition of the field contents.

.. report:: Mapping.PicardDuplicatesMetrics
   :render: xls-table

   Duplicates summary

.. report:: Mapping.PicardDuplicatesHistogram
   :render: line-plot
   :as-lines:

   Duplicates distribution
