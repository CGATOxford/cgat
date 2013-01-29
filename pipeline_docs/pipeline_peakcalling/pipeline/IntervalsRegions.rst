=======
Regions
=======

This page summarizes the :term:`regions` that have been called with the
various region callers. 

``Regions`` are broader areas of enrichment in the genome, arising from more diffuse signals, such as chromatin marks.


Number of intervals
===================

The following figure shows the number of intervals per set:

.. report:: Intervals.RegionsSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: nintervals

   Number of intervals in all runs

The following table presents the number of intervals and 
the average interval width per set:

.. report:: Intervals.RegionsSummaryTable
   :render: table

   Number of intervals and lengths of intervals in
   all runs

Interval sizes
===============

The following figure shows the average interval length per set:

.. report:: Intervals.RegionsSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: avg(length)

   Average interval length

The following plot shows the distribution of the interval size for each set.

.. report:: Intervals.RegionsLengths
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths

Distribution of peak scores
===========================

Average values
--------------

The following plot shows the distribution of the average read coverage within intervals.

.. report:: Intervals.RegionsAverageValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the average read depth within intervals.

Values at peak
--------------

The following plot shows the distribution of the maximum read coverage within intervals.

.. report:: Intervals.RegionsPeakValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the maximum read depth within intervals.

Peak location
=============

The following plot shows the distribution of the peak location within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.RegionsPeakLocation
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :as-lines:
  
   Distance of peak towards start/end of interval normalized
   by the size of the interval.

The following plot shows the distribution of the peak distance within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.RegionsPeakDistance
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
  
   Distance of peak towards start/end of interval


.. Summary of reads under peaks
.. ============================

.. The following tables show the number of reads for each track that fall under peaks in all tracks

.. .. report:: ReadsUnderPeaks.ReadCountSummary
..    :render: matrix
..    :transform-matrix: correspondence-analysis

..    Total number of reads from each track that fall under peaks


.. .. report:: ReadsUnderPeaks.NormalisedTable
..    :render: table

..    Table showing the normalized number of reads falling under peaks for each track


