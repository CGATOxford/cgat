================
Mapping overview
================

This section gives an overview over the mapping process. 

Tophat results
++++++++++++++

Mapping reads with tophat_ across splice-junctions is the first step in the RNASeq 
analysis pipeline. Tophat now implements mapping of reads to a reference transcriptome,
carrying over unmapped reads to genomic mapping.

The following table present an overview of tophat results for 
each :term:`track`.

.. report:: Mapping.TophatSummary
   :render: table

   Tophat results

.. note:: 

   Unmapped reads are not carried over into subsequence analyses, hence the number
   of input and mapped reads will be the same. Only the latter is shown.

Filtering
+++++++++

The RNASeq pipeline offers two filtering options:

1. Remove reads that map better to a reference transcriptome. These are
   often mismapped reads in which the splice-site detection failed.
   Although Tophat implements this step, the option for filtering has been
   left in the pipeline.

2. [Optional] Remove reads that are non-unique. By default these are left in.

+------------------------------+--------------------------------------------------+
|*Column*                      |*Content*                                         |
+------------------------------+--------------------------------------------------+
|input                         |genomic alignments input to the filtering stage   |
+------------------------------+--------------------------------------------------+
|output                        |genomic alignments output after the filtering     |
+------------------------------+--------------------------------------------------+
|removed_mismapped             |genomic alignments removed because they are likely|
|                              |mismapped alignments                              |
+------------------------------+--------------------------------------------------+
|removed_contigs               |genomic alignmentsn removed because aligned to    |
|                              |unwanted contigs (e.g. chrM).                     |
+------------------------------+--------------------------------------------------+
|removed_nonunique_alignments  |genomic alignments that have been removed because |
|                              |they are non-unique                               |
+------------------------------+--------------------------------------------------+
   
.. report:: Mapping.FilteringSummary
   :render: table
   :slices: input,output,removed_mismapped,removed_nonunique_alignments,removed_contigs

   Filtering summary

.. report:: Mapping.FilteringSummary
   :render: interleaved-bar-plot
   :slices: input,output,removed_mismapped,removed_nonunique_alignments,removed_contigs

   Filtering summary

Alignments
++++++++++

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.BamSummary
   :render: table
   :tracks: r(.accepted)
   :slices: mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :tracks: r(.accepted)
   :slices: mapped,reverse,rna,duplicates

   Mapping summary

.. report:: Mapping.MappingFlagsMismatches
   :tracks: r(.accepted)
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of alignments per number of mismatches in alignment.

Reads
+++++

The following table 

.. report:: Mapping.BamSummary
   :render: table
   :tracks: r(.accepted)
   :slices: reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :tracks: r(.accepted)
   :slices: reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping summary

.. report:: Mapping.MappingFlagsHits
   :tracks: r(.accepted)
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of reads per number of alignments (hits) per read.

