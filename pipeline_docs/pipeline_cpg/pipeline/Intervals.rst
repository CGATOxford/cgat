=================
Binding Intervals
=================

Number and Length of Binding Intervals
======================================

The following table presents the number of intervals and 
the average interval length per dataset:

.. report:: Intervals.IntervalsSummary
   :render: table

   Number of intervals and lengths of intervals in all runs

The following plot shows the distribution of interval length for each set.

.. report:: Intervals.IntervalLengths
   :render: box-plot

   Distribution of interval lengths


.. report:: Intervals.IntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths


Binding Interval Coverage
=========================

Average Coverage
----------------

The following plot shows the distribution of average interval coverage for each set.
The average coverage is the average number of reads covering the bins that constitute the interval.

.. report:: Intervals.IntervalAverageValues
   :render: box-plot
   :logscale: y

   Boxplot of the average number of reads within an interval

.. report:: Intervals.IntervalAverageValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,50
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the average number of reads
   matching to bins within an interval.

Maximum Coverage
----------------

The following plot shows the maximum interval coverage for each set.
The maximum coverage is the maximum number of reads falling into the
bins that constitute an interval. The interval peak is the position with the maximum
number of reads.

.. report:: Intervals.IntervalPeakValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.


Peak location
=============

The following plot shows the distribution of the peak location within
an interval for each set, i.e. is it on the sides or the middle
of an interval. Note that this counting does not take into account
strandedness.

.. report:: Intervals.PeakLocation
   :render: line-plot
   :transform: histogram
   :groupby: all
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
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
  
   Distance of peak towards start/end of interval


Correlation of Interval Length and Coverage
===========================================

The following table tests if there is a correlation 
between interval length and average coverage within each set.

.. report:: Intervals.IntervalLengthVsAverageValue
   :render: table
   :transform: correlation
 
   Scatter plots showing the correlation between 
   length and average coverage

The following table tests if there is a correlation 
between interval length and the maximum coverage within each set.

.. report:: Intervals.IntervalLengthVsPeakValue
   :render: table
   :transform: correlation
 
   Scatter plots showing the correlation between 
   length and maximum coverage

