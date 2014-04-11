=====
Peaks
=====

This page summarizes the :term:`peaks` that have been called with the
various peak callers. 

Peaks here refer to intervals arising from the detection of point binding factors, such as classical transcription factors.

Number of intervals
===================

The following figure shows the number of intervals per set:

.. report:: Intervals.PeaksSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: nintervals

   Number of intervals in all runs

The following table presents the number of intervals and 
the average interval width per set:

.. report:: Intervals.PeaksSummaryTable
   :render: table

   Number of intervals and lengths of intervals in
   all runs

Interval sizes
===============

The following figure shows the average interval length per set:

.. report:: Intervals.PeaksSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: avg(length)

   Average interval length

The following plot shows the distribution of the interval size for each set.

.. report:: Intervals.PeaksLengths
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :layout: column-3

   Distribution of interval lengths

Distribution of peak scores
===========================

Average values
--------------

The following plot shows the distribution of the average read coverage within intervals.

.. report:: Intervals.PeaksAverageValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:
   :layout: column-3

   Distribution of the average read depth within intervals.

Values at peak
--------------

The following plot shows the distribution of the maximum read coverage within intervals.

.. report:: Intervals.PeaksPeakValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:
   :layout: column-3

   Distribution of the maximum read depth within intervals.

Peak location
=============

The following plot shows the distribution of the peak location within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.PeaksPeakLocation
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :as-lines:
   :layout: column-3
  
   Distance of peak towards start/end of interval normalized
   by the size of the interval.

The following plot shows the distribution of the peak distance within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.PeaksPeakDistance
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :layout: column-3
  
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


