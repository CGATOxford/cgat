=======================
MACS Diagnostics
=======================

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


