===================================
MACS Interval Parameter Correlation
===================================

This page examines the pairwise relationships between MACS intervals parameters:

1. Length
2. Average coverage
3. Peak coverage
4. Fold change
5. CpG density
6. GC content
7. Distance to closest TSS

Scatter plots
-------------

.. report:: Intervals.IntervalLengthVsAverageValue
   :render: scatter-plot
   :groupby: track
   :xrange: 0,10000
   :yrange: 0,500
   :width: 400
   :layout: column-2
  
   Length vs Average Coverage

.. report:: Intervals.IntervalLengthVsPeakValue
   :render: scatter-plot
   :groupby: track
   :xrange: 0,10000
   :yrange: 0,1000
   :width: 400
   :layout: column-2
  
   Length vs Peak Coverage

.. report:: Intervals.IntervalLengthVsFoldChange
   :render: scatter-plot
   :groupby: track
   :xrange: 0,10000
   :yrange: 0,200
   :width: 400
   :layout: column-2
  
   Length vs Fold Change

.. report:: Intervals.IntervalAvgValVsPeakVal
   :render: scatter-plot
   :groupby: track
   :xrange: 0,500
   :yrange: 0,1000
   :width: 400
   :layout: column-2
  
   Average coverage vs peak coverage

.. report:: Intervals.IntervalAvgValVsFoldChange
   :render: scatter-plot
   :groupby: track
   :xrange: 0,500
   :yrange: 0,200
   :width: 400
   :layout: column-2
  
   Average coverage vs fold change

.. report:: Intervals.IntervalPeakValVsFoldChange
   :render: scatter-plot
   :groupby: track
   :xrange: 0,500
   :yrange: 0,200
   :width: 400
   :layout: column-2
  
   Peak coverage vs fold change


