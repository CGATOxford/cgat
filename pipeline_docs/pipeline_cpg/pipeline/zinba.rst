==================
ZINBA Peak Calling
==================

ZINBA Summary
=============

.. report:: zinba.zinbaSummary
   :render: table

   Table of ZINBA summary statistics (paired input sample)

Overlap with MACS Intervals
===========================

.. report:: zinba.OverlapZinbaMacs
   :render: table

   Overlap of MACS and ZINBA intervals

ZINBA Interval Length
=====================

The following plot shows the distribution of interval length for each set.


.. report:: zinba.zinbaIntervalLengths
   :render: box-plot

   Distribution of interval length

.. report:: zinba.zinbaIntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval length


ZINBA Peak Value
=======================

The following plot shows the distribution of maximum (peak) coverage.

.. report:: zinba.zinbaMaxVal
   :render: box-plot
   :logscale: y

   Boxplot of maximum coverage

.. report:: zinba.zinbaMaxVal
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,500000,10
   :xrange: 0,200
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of maximum coverage


