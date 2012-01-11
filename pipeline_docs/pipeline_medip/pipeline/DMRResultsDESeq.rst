=============
DESeq Results
=============

Differential methylation
========================

Summary
-------

The following table summarizes the results of the test for
differential methylation using :term:`DESeq`.

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

.. report:: DifferentialMethylation.TrackerDMRSummaryDESeq
   :render: table

   Table of DESeq differential expression results

.. report:: DifferentialMethylation.TrackerDMRSummaryDESeq
   :render: interleaved-bar-plot
   :force:
   :slices: significant

   Number of significant tests

.. Pairwise plots
.. --------------

.. The following scatter plots show the fold change versus the expression
.. level. Differences called significant are red, all others are black.

.. .. report:: DifferentialMethylation.TrackerDEPairwiseDESeq
..    :render: user
   
..    Differential expression results

.. DESeq Summary
.. =============

.. .. report:: DifferentialMethylation.TrackerDESummaryPlotsDESeq
..    :render: user
   
..    Summary of differential expression analysis

.. DESeq Fitting overview
.. ======================

.. The following section displays summary stats on the model
.. fit.

.. DESeq attempts to estimate the variance of expression levels
.. for a gene from the mean expression level. As the number
.. of replicates is usually low, the observed variance can vary
.. widely around the estimate. The plots in this section examine
.. the fit:

.. Estimated per-gene variance versus model based estimated of variance.
..    Show are the estimated variance of a gene estimated from this gene alone
..    versus the estimated derived from the local fit for this gene through the base
..    level variances of all genes.

.. Residual ECDF plot
..    cumulative distribution of residuals. These should follow
..    a straight line. The plot is stratified by :term:`base level`.
..    The expression level increases from red to blue.

.. .. report:: DifferentialMethylation.TrackerDESeqFit
..    :render: user
   
..    Fitting overview

.. Volcano plots
.. =============

.. .. report:: DifferentialMethylation.VolcanoPlotDESeq
..    :render: scatter-rainbow-plot
..    :groupby: track
..    :layout: column-3
..    :width: 300

..    Volcano plots (log fold change against -log10 pvalue) 
..    for all pairwise comparisons

.. DESeq Glossary
.. ==============

.. .. glossary::

..    base level
..       normalized expression level. The expression level (usually same tag counts)
..       normalized by the library size to make expression level measurements comparable
..       across experiments.

