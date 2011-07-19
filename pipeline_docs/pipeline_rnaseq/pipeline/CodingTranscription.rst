=======================
Protein coding gene set
=======================

This section summarizes expression of protein coding genes.

Gene set build
==============

The following table contains the filtering results from the gene build.
Note that a transcript might fail several filters.

.. report:: Genemodels.BuildSummary
   :render: table
   
   Number of transcripts failing a certain filter.

Gene set content
----------------

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
   :tracks: refcoding
   :slices: genes,transcripts,single_exon_transcripts,contigs,strands
   :force:

   Summary of gene sets - features

Highest expressed genes
=======================

The following table lists the highest expressed genes in each track
according to :term:`cufflinks`.

.. report:: Expression.ExpressionHighestExpressedGenesDetailed
   :render: table
   :force:
   :groupby: track
   :tracks: refcoding_cuffdiff_isoform

   The ten highest expressed transcripts.

Transcript classification
===========================

The following table lists the number of transcripts with respect to their genomic context.

.. report:: Genemodels.TranscriptClassCounts
   :render: table
   :tracks: r(refcoding)
   :groupby: all

   Genomic classification of lincrna

Differentially expressed genes
==============================

The following table lists the results of the differential expression
comparisons:

.. report:: DifferentialExpression.TrackerDESummaryDESeq
   :render: table
   :tracks: r(refcoding)

   Table of DESeq differential expression results.

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: table
   :tracks: r(refcoding)

   Table of cuffdiff differential expression results.





