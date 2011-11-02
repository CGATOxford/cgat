============
Peak calling
============

Normalization
=============

Peak calling results
=====================

.. report:: PeakCalling.MacsSummary
   :render: table
   :force:

   Summary of MACS results

Filtering
=========

.. report:: PeakCalling.MacsFiltering
   :render: table
   :force:

   Table showing the number of peaks called at a given FDR

Diagnostics
===========

The following plots show the proportion of peaks that are found
if only a proportion of reads are used for peak calling.

.. report:: PeakCalling.MacsDiagnostics                                                                                                                                                                                                      
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


Summary of reads under peaks
============================

The following tables show the number of reads for each track that fall under peaks in all tracks

.. report:: ReadsUnderPeaks.ReadCountSummary
   :render: matrix
   :transform-matrix: correspondence-analysis

   Total number of reads from each track that fall under peaks


.. report:: ReadsUnderPeaks.NormalisedTable
   :render: table

   Table showing the normalized number of reads falling under peaks for each track





