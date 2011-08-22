=======
Binding
=======

Fold change
-----------

The following plot shows the distribution of fold change
comparing a foreground set with the corresponding background
(:term:`Unstim`).

.. report:: Intervals.FoldChange
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Distribution of fold change of maximum reads at peak

.. report:: Intervals.FoldChangeCounts
   :render: interleaved-bar-plot

   Number of interval with at least a 2-fold change compared
   to the unstimulated data.


