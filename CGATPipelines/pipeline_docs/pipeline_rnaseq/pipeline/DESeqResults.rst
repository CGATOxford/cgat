=============
DESeq Results
=============

Differential expression
=============================

Summary
-------

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
|twofold                                 |Number of significant tests, where the  |
|                                        |expression level differs by at least a  |
|                                        |factor of 2.                            |
+----------------------------------------+----------------------------------------+

.. report:: DifferentialExpression.TrackerDESummaryDESeq
   :render: table

   Table of DESeq differential expression results

.. report:: DifferentialExpression.TrackerDESummaryDESeq
   :render: interleaved-bar-plot
   :force:
   :slices: significant

   Number of significant tests

Pairwise plots
--------------

The following scatter plots show the fold change versus the expression
level. Differences called significant are red, all others are black.

.. report:: DifferentialExpression.TrackerDEPairwiseDESeq
   :render: user
   
   Differential expression results

DESeq Summary
=============

.. report:: DifferentialExpression.TrackerDESummaryPlotsDESeq
   :render: user
   
   Summary of differential expression analysis

Volcano plots
=============

.. report:: DifferentialExpression.VolcanoPlotDESeq
   :render: scatter-rainbow-plot
   :groupby: track
   :layout: column-3
   :width: 300

   Volcano plots (log fold change against -log10 pvalue) 
   for all pairwise comparisons

DESeq status 
============

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: export/deseq/*.png
   :width: 200
   :layout: column-3
   
   Various summary plots created by DESeq.

DESeq Glossary
==============

.. glossary::

   base level
      normalized expression level. The expression level (usually same tag counts)
      normalized by the library size to make expression level measurements comparable
      across experiments.

