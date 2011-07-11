====================
Intervals Comparison
====================

Unique Intervals Per Dataset
============================

The following table presents the number of intervals that were uniue to a particular dataset:

.. report:: compareIntervals.UniqueIntervals
   :render: table

   Intervals unique to an individual dataset

Length
------

The following plot shows the distribution of interval length for each set.

.. report:: compareIntervals.UniqueIntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths

Average Coverage
----------------

The following plot shows the distribution of average interval coverage for each set.
The average coverage is the average number of reads covering the bins that constitute the interval.

.. report:: compareIntervals.UniqueIntervalAverageValues
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

.. report:: compareIntervals.UniqueIntervalPeakValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.

Fold Change
-----------

The following plot shows the fold change over control (input) for each set.

.. report:: compareIntervals.UniqueIntervalFoldChange
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of fold enrichment for interval compared to control.


Coverage of Unique Intervals in Different Datasets
--------------------------------------------------

The following chart plots the distribution of interval coverage in datasets 
where binding intervals were not detected.

.. report:: compareIntervals.UniqueIntervalCoverage
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval coverage



Shared Intervals
================

The following table presents the number of intervals in each dataset, 
which overlap with one or more intervals in all other datasets:

.. report:: compareIntervals.SharedIntervals
   :render: table

   Intervals shared between all datasets

Length
------

The following plot shows the distribution of interval length for each set.

.. report:: compareIntervals.sharedIntervalLengths
   :render: line-plot
   :transform: histogram
   :groupby: all
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:

   Distribution of interval lengths

Average Coverage
----------------

The following plot shows the distribution of average interval coverage for each set.
The average coverage is the average number of reads covering the bins that constitute the interval.

.. report:: compareIntervals.SharedIntervalAverageValues
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

.. report:: compareIntervals.SharedIntervalPeakValues
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of the number of reads at the peak within an interval.
   The distribution list the proportion of intervals of a certain peak
   value or more.

Fold Change
-----------

The following plot shows the fold change over control (input) for each set.

.. report:: compareIntervals.SharedIntervalFoldChange
   :render: line-plot
   :transform: histogram
   :groupby: all
   :tf-range: 0,100
   :tf-aggregate: normalized-total,reverse-cumulative
   :as-lines:

   Distribution of fold enrichment for interval compared to control.





