******************
Distances
******************

This part of the pipeline examines the distances of transcript
model to the reference gene set.

Closest distance
================

The following two plots show the distribution of the distance
between transcipt models and the closest gene in the reference
gene set. Transcript models overlapping genes are ignored.

.. report:: Trackers.ClosestDistance
   :transform: histogram
   :render: line-plot
   :logscale: x
   :tf-aggregate: normalized-total
   :tf-bins: log-100
   :as-lines:

   Histogram of distances. The values have been normalized.

.. report:: Trackers.ClosestDistance
   :transform: histogram
   :render: line-plot
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: log-100
   :as-lines:

   Histogram of distances. The values have been normalized
   and the cumulative histogram is shown.
   


