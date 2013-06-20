================
LincRNA gene set
================

The lincrna gene set is constructed by filtering all transcripts
predicted by :term:`cufflinks`. 

Gene set build
==============

The following table contains the filtering results from the gene build.
Note that a transcript might fail several filters.

.. report:: Genemodels.BuildSummary
   :render: table
   :tracks: lincrna
   
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
   :tracks: lincrna
   :slices: genes,transcripts,single_exon_transcripts,contigs,strands
   :force:

   Summary of gene sets - features

Highest expressed lincRNAs
==========================

The following table lists the highest expressed :term:`lincrna` in each track
according to :term:`cufflinks`.

.. report:: Expression.ExpressionHighestExpressedGenesDetailed
   :render: table
   :force:
   :groupby: track
   :tracks: lincrna_cuffdiff_isoform

   The ten highest expressed lincrna transcripts.

LincRNA classes
===============

The following table lists the number of transcripts with respect to their genomic context.

.. report:: Genemodels.TranscriptClassCounts
   :render: table
   :tracks: r(lincrna)
   :groupby: all

   Genomic classification of lincrna

Differentially expressed lincRNAs
=================================

The following table lists the results of the differential expression
comparisons:

.. report:: DifferentialExpression.TrackerDESummaryDESeq
   :render: table
   :tracks: r(lincrna)

   Table of DESeq differential expression results for :term:`lincrna` only

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: table
   :tracks: r(lincrna)

   Table of cuffdiff differential expression results for :term:`lincrna` only.





