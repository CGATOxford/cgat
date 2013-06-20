==================
SICER Peak Calling
==================

SICER Summary
=============

.. report:: sicer.SicerSummary
   :render: table

   Table of SICER summary statistics (paired input sample)

Overlap with MACS Intervals
===========================

.. report:: sicer.OverlapMacs
   :render: table

   Overlap of MACS and SICER intervals

SICER Interval Length
=====================

The following plot shows the distribution of interval length for each set.

.. report:: sicer.SicerIntervalLengths
   :render: box-plot

   Distribution of interval lengths


.. report:: sicer.SicerIntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths


SICER Fold Change
=================

The following plot shows the distribution of interval fold change between input and CAP-seq.

.. report:: sicer.SicerFoldChange
   :render: box-plot
   :logscale: y

   Boxplot of the fold change

.. report:: sicer.SicerFoldChange
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of fold change


