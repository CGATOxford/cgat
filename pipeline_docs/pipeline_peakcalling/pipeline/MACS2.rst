=====
MACS2
=====

This page summarizes MACS2 results.

Normalization
=============

Peak calling results
=====================

.. report:: PeakCalling.Macs2Summary
   :render: table
   :force:

   Summary of MACS results

Filtering
=========

.. report:: PeakCalling.Macs2Filtering
   :render: table
   :force:

   Table showing the number of peaks called at a given FDR

Diagnostics
===========

The following plots show the proportion of peaks that are found
if only a proportion of reads are used for peak calling.

.. report:: PeakCalling.Macs2Diagnostics                                                                                                                                                                                                      
   :render: line-plot                                                                                                                                                                                                                        
   :transform: filter                                                                                                                                                                                                                        
   :tf-fields: proportion of reads,proportion of peaks
   :groupby: track                                                                                                                                                                                                                           
   :as-lines:                                                                                                                                                                                                                                
   :width: 300
   :layout: column-3
   :yrange: 0,110

   Proportion of peaks that are found by MACS if only a proportion of
   peaks are used for peak calling. The different lines correspond to
   peaks of different heights.





