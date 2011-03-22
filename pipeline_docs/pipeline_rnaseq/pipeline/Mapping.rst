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
|:term:`track`.genome.bam               |Alignments of reads mapped on the chromosome      |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.trans.bam                |Alignments of reads reads mapped to transcripts   |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.mismapped.bam            |Alignments flagged as :term:`mismapped`           |
+---------------------------------------+--------------------------------------------------+


Alignments
----------

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: total,mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: total,mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per number of mismatches in alignment.

Reads
-----

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

Tophat results
==============

The following table present an overview of tophat
results for each :term:`track`.

.. report:: Mapping.TophatSummary
   :render: table

   Tophat results

Context results
===============

The following table lists the genomic context that reads map to. Counts are in terms of alignments.
Note that some of these contexts can be overlapping, thus some alignments might be counted several
times. Also, an alignment is assigned to the genomic context that it overlaps by at least 50%. Thus some
alignments spanning several contexts might be dropped.

.. report:: Mapping.MappingContext
   :render: table
   :force:

   Number of alignments that align in a certain genomic context

Ribosomal RNA
-------------

Ribosomal RNA is one of the most abundant transcripts in a cell and dominates RNASeq samples
until it is removed. The following plots and tables examine the number of alignments to
repetitive RNA. Repetetive RNA annotation is taken from the UCSC repeatmasker tracks.

.. report:: Mapping.MappingContext
   :render: table
   :slices: mapped,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA

   Number of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track

.. report:: Mapping.MappingContext
   :render: pie-plot
   :pie-first-is-total: notRNA
   :groupby: track
   :slices: mapped,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA
   :layout: column-3
   :width: 200

   Proportion of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track


Protein coding expression
-------------------------

The following plots list the number of alignments to protein coding and (protein coding) 
pseudogene exons. The annotations are taken from the ENSEMBL gene set.

.. report:: Mapping.MappingContext
   :render: pie-plot
   :pie-first-is-total: genomic
   :groupby: track
   :slices: mapped,protein_coding,pseudogene
   :layout: column-3
   :width: 200

   Proportion of alignments that align to protein coding genes or pseudo genes.






