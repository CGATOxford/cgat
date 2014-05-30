=====================================
Exported interval sets
====================================

Filtered interval sets
======================

Intervals are output after filtering. Filtering are removing intervals
that

* overlap regions with a large background signal as determined from
  input data

* overlap with each other. The final interval set will contain only
  non-overlapping intervals.

The table below displays the filtering summary for each data set.

.. report:: PeakCalling.FilteringSummary
   :render: stacked-bar-plot
   :slices: output,removed_background,removed_merged
   :split-at: 20
   :layout: column-2

   Filtering overview
	    
.. report:: PeakCalling.FilteringSummary
   :render: table

   Table with export summary

                   
