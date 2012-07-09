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

FastQC
======

The following list presents links to the results of the :term:`fastqc` tool.
The fastqc tool is run over aligned reads only.

.. report:: Mapping.FastQCReport
   :render: user

   fastqc results


Picard
======

The following table present an overview of the alignments in the 
BAM files for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.PicardSummary
   :tracks: r(.accepted)
   :render: table

   Alignments summary

.. report:: Mapping.PicardSummary
   :tracks: r(.accepted)
   :render: interleaved-bar-plot
   :slices: PCT_PF_READS_ALIGNED,STRAND_BALANCE

   Percentage quantities

.. report:: Mapping.PicardSummary
   :tracks: r(.accepted)
   :render: interleaved-bar-plot
   :slices: PF_READS_ALIGNED,PF_HQ_ALIGNED_READS

   Percentage quantities

.. report:: Mapping.AlignmentQualityByCycle
   :tracks: r(.accepted)
   :render: line-plot
   :as-lines:
   :yrange: 0,

   mean quality score by cycle

.. report:: Mapping.AlignmentQualityDistribution
   :tracks: r(.accepted)
   :render: line-plot
   :as-lines:
   :yrange: 0,

   quality score distribution

