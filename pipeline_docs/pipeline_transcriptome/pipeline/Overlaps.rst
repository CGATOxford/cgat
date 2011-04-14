.. _OverlapBetweenSets:

*****************************
Overlap between sets
*****************************

This page presents the overlap between sets.

Number of transcript models present in all sets
===============================================

.. report:: Trackers.SharedGenes                                                                                                                                                                                                             
   :render: table
   :groupby: track

   Transcript models present in all sets. Computation is based on the set
   ``merged``, thus there is a possibility that in a pairwise comparison,
   two transcript models might not overlap.

Plots of pairwise overlap between sets
======================================

The following table :ref:`PipeFigureOverlapBetweenSetsPercent` lists 
the percent of genes in a track that overlap with genes in another track.
Each entry is computed as the ratio of ``number of overlapping genes between
row and column`` and the ``number of genes in the set row``. Hence, the
matrix is not symmetric.

.. _PipeFigureOverlapBetweenSetsPercent:

.. report:: Trackers.OverlapsGenesPercent
   :render: matrix-plot

   Percent overlap between sets. The matrix plot displays for each set in a row 
   the percentage of that transcripts models that overlap with transcript models
   in the set in the column. 

Note that in the reference gene set (often called *ensembl*) only coding 
exons are included in the analysis. Thus in the slice *known*, not necessarily all 
transcript models will overlap with the reference set as pure UTR
transcript models are excluded.

Tables of overlap between sets
==============================

The following tables list the percentage of genes in a track that overlap with genes in another track.
Each entry is computed as the ratio of ``number of overlapping row genes between
row and column`` and the ``number of genes in the set row``. Hence, the
table is not symmetric.

.. _PipeTableOverlapBetweenSetsPercent:

.. report:: Trackers.OverlapsGenesPercent
   :render: matrix

   Percent overlap between sets counted by genes

The following tables list the number of genes in a track that overlap with genes in another track.
Note that the number of overlapping genes between two sets will differ unless they contain the
exact same gene models. For example, if a gene A in set 1 is split into genes A1 and A2 in set2,
the number of overlapping genes will be 1 and 2 for set 1 and set 2, respectively. Hence, the table
is not symmetric.

.. _PipeTableOverlapBetweenSetsCounts:

.. report:: Trackers.OverlapsGenesCounts
   :render: matrix

   Overlap between gene sets as counts

The following tables list the number of bases in a track that overlap with genes in another track.
Each entry is computed as the ratio of ``number of overlapping row genes between
row and column`` and the ``number of genes in the set row``. Hence, the
table is not symmetric.

.. report:: Trackers.OverlapsBasesPercent
   :render: matrix

   Percent base overlap between sets

Annotator significance of overlap
=================================

Significance of overlap between transcript models.

.. report:: Annotator.AnnotatorMatrixSets
   :render: matrix-plot
   :transform-matrix: correspondence-analysis,transpose
   :palette: RdBu
   :zrange: -1000,1000

   Statistical significant enrichment (in percent) computed by annotator between novel
   transcript models. Plotted are fold enrichment (positive values) and 
   fold depletion (negative values). The maxmimum enrichment/depletion shown is 
   10-fold. The workspace in this analysis is the full genome.

.. report:: Annotator.AnnotatorFullTableSets
   :render: table

   Full table with annotator results (including non-signficant ones)
