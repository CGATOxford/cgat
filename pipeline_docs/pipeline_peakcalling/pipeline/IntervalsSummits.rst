=======
Summits
=======

This page summarizes the :term:`region` that have been called with the
various peak callers. 

``Regions`` are larger areas in the genome with higher read density
than expected as opposed to ``summits``, which are smaller summits 
within regions.

Number of intervals
===================

The following figure shows the number of intervals per set:

.. report:: Intervals.SummitsSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: nintervals

   Number of intervals in all runs

The following table presents the number of intervals and 
the average interval width per set:

.. report:: Intervals.SummitsSummary
   :render: table

   Number of intervals and lengths of intervals in
   all runs

Interval sizes
===============

The following figure shows the average interval length per set:

.. report:: Intervals.SummitsSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: avg(length)

   Average interval length

The following plot shows the distribution of the interval size for each set.

.. report:: Intervals.SummitsLengths
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

.. report:: Intervals.SummitsAverageValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the average read depth within intervals.

Values at peak
--------------

The following plot shows the distribution of the maximum read coverage within intervals.

.. report:: Intervals.SummitsPeakValues
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

.. report:: Intervals.SummitsPeakLocation
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

.. report:: Intervals.SummitsPeakDistance
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
  
   Distance of peak towards start/end of interval

