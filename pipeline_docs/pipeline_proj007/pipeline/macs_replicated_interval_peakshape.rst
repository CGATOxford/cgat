=====================================================
Heatmap of Liver / Testes Specific CAPseq intervals
=====================================================

This figure shows the normalised read density in a 5kb window centered on CAPseq intervals. 

The figures are scaled from 0 to 100 and not from 0 to max in order 
to increase the visibility of lower coverage tracks.

The syntax for track names is ``<bam-files>.<bed-file>.peakshape``

.. report:: peakshape.PeakShapeTrackerCentre
   :render: matrixNP-plot
   :slices: interval_score
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,0.5
   
   Heatmaps of merged liver and testes intervals sorted by tissue specificity then interval width (centred on peak middle)
   
.. report:: peakshape.PeakShapeTracker
   :render: matrixNP-plot
   :slices: interval_score
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,0.5
   
   Heatmaps of merged liver and testes intervals sorted by tissue specificity then interval width (centred on highest peak)
   
.. report:: peakshape.PeakShapeTrackerCentreNoScale
   :render: matrixNP-plot
   :slices: interval_score
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,0.5
   
   Heatmaps of merged liver and testes intervals sorted by tissue specificity then interval width (centred on peak middle) with no subsampling
   
   																												
