.. _GOAnnotatorsTerritories:

*************
GO annotation
*************

This part of the pipeline uses Gerton Lunter's ``annotator``
program to check for enrichment or depletion of GO annotations.
GO annotations are associated with genomic segments through
genes. The genomic location of each gene from the first to
the last exon defines its ``territory``. The territory
is extended by 30kb in either direction at most. Overlapping
territories are resolved by assigning the boundary to me the
mid-point between the terminal exons of the two neighbouring
genes.

The output uses a false discovery rate threshold of q=0.05.

GO territories
==============

Overview
--------

.. report:: Trackers.AnnotatorSummaryGOSlimTerritories
   :render: table
   :groupby: all

   Counts of over/under represented GOSlim categories.

.. report:: Trackers.AnnotatorMatrixGOSlimTerritories
   :render: matrix-plot
   :zrange: -100,100
   :transform-matrix: correspondence-analysis,transpose
   :palette: RdBu

   Visual overview of over- or under-represented GOSlim categories
   in blue and red, respectively. The scale is fold enrichment/
   depletion in percent over the expected. Shown are weakly
   enriched/depleted categories (<100%).
   Categories and samples without enrichment or depletion are
   not shown.

.. report:: Trackers.AnnotatorMatrixGOSlimTerritories
   :render: matrix-plot
   :zrange: -400,400
   :transform-matrix: correspondence-analysis,transpose
   :palette: RdBu

   Visual overview of over- or under-represented GOSlim categories
   in blue and red, respectively. The scale is fold enrichment/
   depletion in percent over the expected. Shown are strongly
   enriched/depleted categories (<400%).
   Categories and samples without enrichment or depletion are
   not shown.

Detailed categories
-------------------

Table with overrepresented GOSlim

.. report:: Trackers.AnnotatorEnrichmentGOSlimTerritories
   :render: grouped-table

   GOSlim annotations of transcripts models. Only enriched categories are shown.

Table with underrepresented GOSlim

.. report:: Trackers.AnnotatorDepletionGOSlimTerritories
   :render: grouped-table

   GOSlim annotations of transcript models. Only depleted categories are shown.

