=========================
Heatmap of CAPseq signal
=========================

This figure shows the read density in various intervals.  The
intervals have been sorted by peak width. The further down in the
image, the narrower the peak. 

The figures are scaled from 0 to 100 and not from 0 to max in order 
to increase the visibility of lower coverage tracks.

The syntax for track names is ``<bam-files>.<bed-file>.peakshape``

.. report:: peakshape.PeakShapeTracker
   :render: matrixNP-plot
   :slices: interval_width
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,100
   
   Heatmaps of merged liver and testes intervals sorted by interval width
   
.. report:: peakshape.PeakShapeTracker
   :render: matrixNP-plot
   :slices: interval_score
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,100
   
   Heatmaps of merged liver and testes intervals sorted by tissue specificity then interval width
   																													
.. report:: peakshape.PeakShapeTracker
   :render: matrixNP-plot
   :slices: peak_width
   :layout: column-2
   :palette: Blues
   :nolabel-cols:
   :max-rows: 9
   :nolabel-rows:
   :max-cols: 0
   :zrange: 0,100
   
   Heatmaps of merged liver and testes intervals sorted by peak width
   
   																													



