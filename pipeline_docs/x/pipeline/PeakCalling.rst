============
Peak calling
============


MACS Peak calling results
=========================

.. report:: PeakCalling.MacsSummary
   :render: table

Diagnostics
===========

The following plots show the proportion of peaks that are found
if only a proportion of reads are used for peak calling.

.. report:: PeakCalling.MacsDiagnostics
   :render: line-plot
   :transform: filter
   :tf-fields: proportion of reads, proportion of peaks
   :groupby: track
   :as-lines:
   :width: 300
   :layout: column-3
   :yrange: 0,110

   Proportion of peaks that are found by MACS if only a proportion of
   reads are used for peak calling. The different lines correspond to
   peaks of different heights.

