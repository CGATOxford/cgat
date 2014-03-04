====================
Alignment statistics
====================

This section brings together various alignment statistics.

.. Bamstats
.. ========

.. The following list presents links to the results of the :term:`bamstats` tool.

.. .. report:: Mapping.BamReport
..    :render: user

..    bamstats results

.. FastQC
.. ======

.. The following list presents links to the results of the :term:`fastqc` tool.
.. The fastqc tool is run over aligned reads only.

.. .. report:: Mapping.FastQCReport
..    :render: user

..    fastqc results


Picard metrics
==============

The following table present an overview of the alignments in the 
BAM files for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.PicardSummary
   :render: table

   Alignments summary

.. report:: Mapping.PicardSummary
   :render: interleaved-bar-plot
   :slices: PCT_PF_READS_ALIGNED,STRAND_BALANCE
   :layout: row
   :split-at: 10

   Percentage quantities

.. report:: Mapping.PicardSummary
   :render: interleaved-bar-plot
   :slices: PF_READS_ALIGNED,PF_HQ_ALIGNED_READS
   :layout: row
   :split-at: 10

   Percentage quantities

The following plots show the distribution of base quality scores in
reads. Note that if the bam-files have been stripped of sequence and
quality information the plots below will contain meaningless values.

.. report:: Mapping.AlignmentQualityByCycle
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :layout: row

   mean quality score by cycle

.. report:: Mapping.AlignmentQualityDistribution
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :layout: row

   quality score distribution

The following table presents an overview of the insert size metrics
for each :term:`track`.  See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#InsertSizeMetrics>`_
for a definition of the field contents.
(note: only applicable for paired-end data).

.. report:: Mapping.PicardInsertSizeMetrics
   :render: table

   Insert size summary

.. report:: Mapping.PicardInsertSizeHistogram
   :render: line-plot
   :as-lines:

   Insert size distribution

Coverages
=========

The following base coverage summaries are computed from bigwig files.
This report might be empty if no bigwig files have been computed.

.. report:: Mapping.BigwigSummary
   :render: table
   :force:

   Coverage summary computed from bigwig files.

