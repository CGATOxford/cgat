==================
Signal/Noise ratio
==================

This section examines the signal/noise ratio in the data set.

.. _snratiomedian:

S/N Ratio based on median
=========================

The following plots take a subset of genomic windows for estimating
the signal/noise ratio.

To create the plots, the read count in each window is divided by the
median read count across all windows. The idea is to observe a subset
of windows that have read-counts that are much higher than the
median. It is helpful to compare these data with input tracks.

The assumption in these plots is that less than half the windows are
enriched for signal and thus the median reflects the number of reads
expected as background.

The first set of plots display the distribution of log2 fold values
as cumulative histograms:

A pseudo-count of 1 is added to every window.

.. report:: Signal.SignalMedian
   :render: line-plot
   :transform: histogram
   :split-at: 10
   :split-always: input,Input
   :as-lines:
   :tf-aggregate: normalized-total,cumulative
   :layout: column-2
   :width: 200
   
   Distribution of read counts in windows
   normalized by median read count within a sample

The same as violin plots:
  
.. report:: Signal.SignalMedian
   :render: violin-plot
   :tight:

   Violin-plot of read counts in windows
   normalized by median read count within a sample

.. _snratioinput:

S/N Ratio based on input
=========================

The following plots take a subset of windows and divide the 
read-counts in each window with the read-count of the same
window in the appropriate input track.

A pseudo-count of 1 is added to every window and the fold-changes
are normalized by the median read count per window per track and input track.

If no input tracks are available, the fold change is computed 
with the per-track median read-count as above.

.. report:: Signal.SignalInput
   :render: line-plot
   :transform: histogram
   :split-at: 10
   :split-always: input,Input
   :as-lines:
   :tf-aggregate: normalized-total,cumulative
   :layout: column-2
   :width: 200
   
   Distribution of read counts in windows
   normalized by median read count within a sample

The same as violin plots:
  
.. report:: Signal.SignalInput
   :render: violin-plot
   :tight:

   Violin-plot of read counts in windows
   normalized by median read count within a sample



