=========
Intervals
=========

Number of intervals in each set
=================================

The following table presents the number of intervals and 
the average interval width per set:

.. report:: Intervals.IntervalsSummary
   :render: table

   Number of intervals and lengths of intervals in all runs

The following figure shows the number of intervals and
average interval width per set:

.. report:: Intervals.IntervalsSummary
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: nintervals,<length>

   Number of intervals and lengths of intervals in all runs

Distribution of peak scores
===========================

Average values
--------------

The following plot shows the average interval value for each set.
The average value is the average number of reads falling into the
bins that constitute an interval.

.. report:: Intervals.IntervalAverageValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the average number of reads
   matching to bins within an interval.

Values at peak
--------------

The following plot shows the maximum (peak) interval value for each set.
The peak value is the maximum number of reads falling into the
bins that constitute an interval. The peak is the position with the maximum
number of reads.

.. report:: Intervals.IntervalPeakValues
   :render: line-plot
   :transform: histogram
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.

Distribution of interval size
=============================

The following plot shows the distribution of interval size for each set.

.. report:: Intervals.IntervalLengths
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths

Peak location
=============

The following plot shows the distribution of the peak location within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.PeakLocation
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

.. report:: Intervals.PeakDistance
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
  
   Distance of peak towards start/end of interval

Correlation of interval size and peak score
===========================================

The following table tests if there is a correlation 
between interval size and average value within each set.

.. report:: Intervals.IntervalLengthVsAverageValue
   :render: table
   :transform: correlation
 
   Scatter plots showing the correlation between 
   length and average value

The following table tests if there is a correlation 
between interval size and the peak value within each set.

.. report:: Intervals.IntervalLengthVsPeakValue
   :render: table
   :transform: correlation
 
   Scatter plots showing the correlation between 
   length and peak value

