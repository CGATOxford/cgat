===========================
Correlation of peak heights
===========================

The following tables and plots show the correlation of
peak heights between various tracks.

.. report:: Intervals.CorrelationsPeakval
   :render: table
   :transform: select,correlation
   :tf-fields: peakval

   Full table

.. report:: Intervals.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :tf-fields: peakval,coefficient
   :transform-matrix: correspondence-analysis
   :zrange: -1,1

   Correlation coefficient

.. report:: Intervals.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :transform-matrix: correspondence-analysis
   :tf-fields: peakval,pvalue

   P-Value

.. report:: Intervals.CorrelationsPeakval                                                                                                                                                                                                     
   :render: scatter-plot                                                                                                                                                                                                                     
   :transform: combine                                                                                                                                                                                                                       
   :tf-fields: peakval                                                                                                                                                                                                                       
   :groupby: none
   :width: 200
   :layout: column-5

   Scatter plots of pairwise combination of peakval.

