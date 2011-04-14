=======
Binding
=======

Fold change
-----------

The following plot shows the distribution of fold change
comparing a foreground set with the corresponding background
(:term:`Unstim`).

.. report:: ChipseqReport.FoldChange
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :tracks: merged

   Distribution of fold change of maximum reads at peak

.. report:: ChipseqReport.FoldChangeCounts
   :render: interleaved-bar-plot
   :tracks: merged

   Number of interval with at least a 2-fold change compared
   to the unstimulated data.


