.. _Cuffdiff results:

=================
Cuffdiff results
=================

Cuffdiff differential expression
=================================

The following table summarizes the results of the test for differential expression.
It lists for every combination of object (gene, isoform, ...), geneset and pair of tracks
the number of tests that passed, failed, etc.

+----------------------------------------+----------------------------------------+
|*Column*                                |*Content*                               |
+----------------------------------------+----------------------------------------+
|significant                             |Number of significant tests             |
+----------------------------------------+----------------------------------------+
|tested                                  |Number of tests in total                |
+----------------------------------------+----------------------------------------+
|status_OK                               |Number of successful tests              |
+----------------------------------------+----------------------------------------+
|status_FAIL                             |Number of failed tests                  |
+----------------------------------------+----------------------------------------+
|status_NOTEST                           |Number of tests ignored                 |
+----------------------------------------+----------------------------------------+
|status_NOCALL                           |Number of successful tests, but in which|
|                                        |at least one gene has a very low        |
|                                        |expression level.                       |
+----------------------------------------+----------------------------------------+
|twofold                                 |Number of significant tests, where the  |
|                                        |expression level differs by at least a  |
|                                        |factor of 2.                            |
+----------------------------------------+----------------------------------------+

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: table
   :force:

   Table of differential expression results

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: interleaved-bar-plot
   :force:
   :slices: significant

   Number of significant tests

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: interleaved-bar-plot
   :force:
   :slices: tested,significant,twofold

   add caption here

The following scatter plots show the fold change versus the expression
level. Differences called significant are red, all others are black.

.. report:: DifferentialExpression.TrackerDEPairwiseCuffdiff
   :render: user
   
   Differential expression results

Cuffdiff Summary
=================

The following plots show p-value versus length. Without length normalization,
longer genes are more likely to be called differentially expressed.

.. report:: DifferentialExpression.TrackerDESummaryPlotsCuffdiff
   :render: user
   
   Summary of differential expression analysis

Volcano plots
=============

.. report:: DifferentialExpression.VolcanoPlotCuffdiff
   :render: scatter-rainbow-plot
   :groupby: track
   :layout: column-3
   :width: 300

   Volcano plots (log fold change against -log10 pvalue) 
   for all pairwise comparisons


 
