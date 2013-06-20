===========================
Sequencing real estate
===========================

Highly expressed genes might take up a lot of the sequencig
real-estate in an experiment. As a consequence, lower expressed transcripts
might not be very well sampled and thereby complicating their analysis.

Genes have been ranked in descending order by the number of reads aligning to them and are
sorted on the x-axis. The y-axis displays the cumulative proportion of reads that
are mapping to these genes.

.. report:: Expression.TrackerRealEstate
   :render: line-plot
   :tracks: r(^(?!.*agg).*$)
   :as-lines:

   Sequencing real estate

Reading from the left, the plot indicates, what proportion of reads
map to a certain number of highest expressed genes. 


