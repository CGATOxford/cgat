===========================
Sequencing real estate
===========================

Highly expressed genes might take up a lot of the sequencig
real-estate in an experiment. As a consequence, lower expressed transcripts
might not be very well sampled and thereby complicating their analysis.

Genes have been ranked in order by the number of reads aligning to them and are
sorted on the x-axis in reverse order with highest expressed genes on
the left of the x-axis and lowest expressed genes on the right of the x-axis.

The y-axis displays the cumulative proportion of reads that are
mapping to these genes. The faster the curve rise, the fewer genes are
responsible for most of the reads.

.. report:: Expression.TrackerRealEstate
   :render: line-plot
   :as-lines:

   Sequencing real estate

Reading from the left, the plot indicates, what proportion of reads
map to a certain number of highest expressed genes. 


