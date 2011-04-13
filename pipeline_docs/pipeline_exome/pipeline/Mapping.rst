=======
Mapping
=======

Mapping results
===============

+---------------------------------------+--------------------------------------------------+
|*Filename*                             |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.bam                      |Alignments after QC and filtering. This is the set|
|                                       |used for subsequent analyses.                     |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.roi.bam                  |Alignments of reads mapped to regions of interest |
+---------------------------------------+--------------------------------------------------+

Reads
-----
.. todo::
   alignment and read tables identical

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.MappingFlagsHits
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of reads per number of alignments (hits) per read.


Alignments
----------

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: total,mapped,reverse,duplicates

   Mapping summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: total,mapped,reverse,duplicates

   Mapping summary

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per number of mismatches in alignment.



Alignment statistics
====================

The following table present an overview of the alignments in the 
BAM files for each :term:`track`. See the 
`Picard metrics <http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics>`_
for a definition of the field contents.
(note: PF=pass filter, reads that pass the vendor's filter criteria).

.. report:: Mapping.AlignmentSummary
   :render: table

   Alignments summary

.. report:: Mapping.AlignmentSummary
   :render: interleaved-bar-plot
   :slices: PCT_PF_READS,PCT_PF_READS_ALIGNED,STRAND_BALANCE

   Percentage quantities

.. report:: Mapping.AlignmentSummary
   :render: interleaved-bar-plot
   :slices: TOTAL_READS,PF_READS,PF_READS_ALIGNED,PF_HQ_ALIGNED_READS

   Percentage quantities

.. report:: Mapping.AlignmentQualityByCycle
   :render: line-plot
   :as-lines:
   :yrange: 0,

   mean quality score by cycle

.. report:: Mapping.AlignmentQualityDistribution
   :render: line-plot
   :as-lines:
   :yrange: 0,

   quality score distribution


