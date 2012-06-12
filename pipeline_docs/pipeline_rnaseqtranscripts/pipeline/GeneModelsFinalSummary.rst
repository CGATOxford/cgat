====================
Summary of gene sets
====================

This page summarizes information about the various genesets built by the
RNASeq pipeline. 

The gene sets of interest are

   * abinitio: the complete ab-initio gene set after removing unique transfrags.
   * lincrna: non-coding transcripts derived from merging the known lincRNA
           transcripts with those derived from the experiment.
   * novel: novel transcripts - transcripts not overlapping any of the coding or non-coding transcripts
           in the :term:`reference` gene set. 
   * refnoncoding: lincRNA and other non-coding transcripts obtained from the the :term:`reference` gene set.
   * refcoding: protein coding transcripts of the the :term:`reference` gene set.
   * reference: the :term:`reference` gene set.

Genes and transcripts
=====================

+------------------------------+--------------------------------------------------+
|*Column*                      |*Content*                                         |
+------------------------------+--------------------------------------------------+
|genes                         |number of genes                                   |
+------------------------------+--------------------------------------------------+
|transcripts                   |number of transcripts                             |
+------------------------------+--------------------------------------------------+
|single_exon_transcripts       |number of single exon transcripts                 |
+------------------------------+--------------------------------------------------+
|features                      |number of features (exons, CDS, ...)              |
+------------------------------+--------------------------------------------------+
|sources                       |number of sources (programs)                      |
+------------------------------+--------------------------------------------------+
|strands                       |number of strands                                 |
+------------------------------+--------------------------------------------------+
|contigs                       |number of chromosomes/contigs                     |
+------------------------------+--------------------------------------------------+

.. report:: Genemodels.GenesetSummary
   :render: table
   :slices: genes,transcripts,single_exon_transcripts,contigs,strands
   :force:

   Summary of gene sets - features

.. report:: Genemodels.GenesetSummary
   :render: interleaved-bar-plot
   :slices: genes,transcripts,single_exon_transcripts
   
   Number of genes/transcripyts

Transcript classes
==================

Transcript classes
------------------

.. report:: Genemodels.TranscriptClassCountsSummaryByClass
   :render: interleaved-bar-plot
   :layout: column-3
   :width: 300
   :transform: select
   :tf-fields: ntranscripts

   Classes of transcripts

.. report:: Genemodels.TranscriptClassCountsSummaryByClass
   :render: table
   :transform: select
   :tf-fields: ntranscripts

   Classes of transcripts

Protein coding classes
----------------------

The following plot lists the various classes of transcripts overlapping protein coding genes.

.. report:: Genemodels.TranscriptClassCounts
   :render: interleaved-bar-plot
   :slices: protein_coding
   :layout: column-3
   :width: 300
   :transform: select
   :tf-fields: ntranscripts

   Overlap with protein coding transcripts

.. report:: Genemodels.TranscriptClassCounts
   :render: table
   :slices: protein_coding
   :transform: select
   :tf-fields: ntranscripts

   Overlap with protein coding transcripts

Exon size summary
=================

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Content*                                         |
+---------------------------------------+--------------------------------------------------+
|exon_count_mean                        |mean number of exons per transcript               |
+---------------------------------------+--------------------------------------------------+
|exon_count_median                      |median number of exons per transcript             |
+---------------------------------------+--------------------------------------------------+
|exon_count_min                         |smallest number of exons per transcript           |
+---------------------------------------+--------------------------------------------------+
|exon_count_max                         |largest number of exons per transcript            |
+---------------------------------------+--------------------------------------------------+
|exon_size_nval                         |number of exons                                   |
+---------------------------------------+--------------------------------------------------+
|exon_size_mean                         |mean exon size                                    |
+---------------------------------------+--------------------------------------------------+
|exon_size_median                       |median exon size                                  |
+---------------------------------------+--------------------------------------------------+
|exon_size_min                          |smallest exon size                                |
+---------------------------------------+--------------------------------------------------+
|exon_size_max                          |largest exon size                                 |
+---------------------------------------+--------------------------------------------------+
 
.. report:: Genemodels.GenesetSummary
   :render: table
   :slices: exon_size_nval,exon_size_mean,exon_size_median,exon_size_min,exon_size_max
   :force:

   Summary of gene sets - exons

.. report:: Genemodels.GenesetSummary
   :render: interleaved-bar-plot
   :slices: exon_size_mean,exon_size_median
   
   Mean/median exon size

.. report:: Genemodels.GenesetSummary
   :render: table
   :slices: exon_count_mean,exon_count_median,exon_count_min,exon_count_max
   :force:

   Summary of gene sets - exons

.. report:: Genemodels.GenesetSummary
   :render: interleaved-bar-plot
   :slices: exon_count_mean,exon_count_median

   Mean/median exon counts

Intron size summary
====================

+------------------------------------------+----------------------------------------------------+
| *Column*                                 |*Content*                                           |
+------------------------------------------+----------------------------------------------------+
| intron_size_nval                         |number of introns                                   |
+------------------------------------------+----------------------------------------------------+
| intron_size_mean                         |mean intron size                                    |
+------------------------------------------+----------------------------------------------------+
| intron_size_median                       |median intron size                                  |
+------------------------------------------+----------------------------------------------------+
| intron_size_min                          |smallest intron size                                |
+------------------------------------------+----------------------------------------------------+
| intron_size_max                          |largest intron size                                 |
+------------------------------------------+----------------------------------------------------+

.. report:: Genemodels.GenesetSummary
   :render: table
   :slices: intron_size_nval,intron_size_mean,intron_size_median,intron_size_min,intron_size_max
   :force:

   Summary of gene sets - introns

