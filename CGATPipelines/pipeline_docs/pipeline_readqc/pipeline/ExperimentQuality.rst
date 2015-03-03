===========================
Per Experiment Read Quality
===========================

This page presents a summary of the distribution of read qualities
aggregated by experiment.

.. report:: ReadQuality.PerExperimentSequenceQuality
   :render: r-ggplot
   :statement: aes(x=Quality, y=value) + geom_histogram(stat="identity") + facet_wrap(~variable) + theme_bw()

   Distribution of read qualities aggregated by experminent.
