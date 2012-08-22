=================
MACS Peak Calling
=================

MACS Summary
============

Summary of the number of peaks (binding intervals) called by MACS.

.. report:: macs.MacsSummary
   :render: table

   Table of MACS summary statistics (paired input sample)


.. report:: macs.MacsSoloSummary
   :render: table

   Table of MACS summary statistics (no control)

MACS Intervals
==================

The following table presents the number of intervals per dataset after filtering with FDR <1%:

.. report:: macs.MacsIntervalsSummary
   :render: table

   MACS interval summary

The following plot indicates the number of MACS intervals that pass a particular fold change threshold.

.. report:: macs.FoldChangeThreshold
   :render: line-plot

   Chart of intervals with fold change exceeding threshold

Background Binding
==================

The following table presents the proportion of reads that overalap binding intervals in each dataset.

.. report:: macs.BackgroundSummary
   :render: table

   Table of the proportion of reads that overlap binding intervals


MACS Diagnostics
================

The following plots show the proportion of peaks that are found
if only a proportion of reads are used for peak calling.

.. report:: macs.MacsDiagnostics
   :render: line-plot
   :transform: filter
   :tf-fields: proportion of reads,proportion of peaks
   :groupby: track
   :as-lines:
   :width: 400
   :layout: column-2
   :yrange: 0,110

   Proportion of peaks that are found by MACS if only a proportion of
   reads are used for peak calling. 

