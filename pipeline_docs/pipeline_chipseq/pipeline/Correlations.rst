************************************
Full correlation matrix of all rates
************************************

Peakval
-------

.. report:: ChipseqReport.CorrelationsPeakval
   :render: table
   :transform: select,correlation
   :tf-fields: peakval

   Full table

.. report:: ChipseqReport.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :tf-fields: peakval,coefficient
   :transform-matrix: correspondence-analysis
   :zrange: -1,1

   Correlation coefficient

.. report:: ChipseqReport.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :transform-matrix: correspondence-analysis
   :tf-fields: peakval,pvalue

   P-Value

.. report:: ChipseqReport.CorrelationsPeakval                                                                                                                                                                                                     
   :render: scatter-plot                                                                                                                                                                                                                     
   :transform: combine                                                                                                                                                                                                                       
   :tf-fields: peakval                                                                                                                                                                                                                       
   :groupby: none
   :width: 200
   :layout: column-5

   Scatter plots of pairwise combination of peakval.

Peakval
-------

.. report:: ChipseqReport.CorrelationsPeakval
   :render: table
   :transform: select,correlation
   :tf-fields: peakval
   :tracks: replicates

   Full table

.. report:: ChipseqReport.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :tf-fields: peakval,coefficient
   :transform-matrix: correspondence-analysis
   :zrange: -1,1
   :tracks: replicates

   Correlation coefficient

.. report:: ChipseqReport.CorrelationsPeakval
   :render: matrix-plot
   :transform: select,correlation,select
   :transform-matrix: correspondence-analysis
   :tf-fields: peakval,pvalue
   :tracks: replicates

   P-Value

.. report:: ChipseqReport.CorrelationsPeakval                                                                                                                                                                                                     
   :render: scatter-plot                                                                                                                                                                                                                     
   :transform: combine                                                                                                                                                                                                                       
   :tf-fields: peakval                                                                                                                                                                                                                       
   :groupby: none
   :width: 200
   :layout: column-5
   :tracks: replicates
   
   Scatter plots of pairwise combination of peakval.

