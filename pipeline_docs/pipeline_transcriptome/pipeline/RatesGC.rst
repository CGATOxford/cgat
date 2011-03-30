*************************
GC content of transcripts
*************************

Statistics
..........

.. report:: Trackers.RatesGC
   :render: table
   :transform: stats

   Statistics of G+C content of transcript models

Plots per slice
...............

The following plots show the distribution of G+C content
in transcript models.

.. report:: Trackers.RatesGC
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,1.1,0.01)
   :as-lines:

   Distribution of G+C content in transcript models
