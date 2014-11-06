=====================
Correlations, heatmaps and clustering
=====================

Global correlation
==========================
The following plots are correlations between treatment groups for all
CpGs, split by CpG island / Non-CpG island.

.. report:: SummaryPlots.correlations
   :render: gallery-plot

   High resolution plots can be downloaded using the links below each plot

Spearman's correlation heatmaps
==========================

The following plots are heatmaps of spearman's rho values for all
possible pairwise comparisons. The samples are ordered by treatment
group and are not clustered. For clustering see the plots in the next
section. Heatmaps are split by CpG island / Non-CpG island.

.. report:: SummaryPlots.heatmaps
   :render: gallery-plot
 
   High resolution plots can be downloaded using the links below each plot


Hierachical clustering
==========================

The following plots are hierarchical clusterings of the samples based
on pairwise spearman's rank correlations. Clusterings are split by are
split by CpG island / Non-CpG island. Distance between samples is
represented by height on the y-axis. The clustering has been
bootstrapped to estimate how strongly the data support the
cluster. The numbers in green are the bootstrap probabilities.

.. report:: SummaryPlots.clustering
   :render: gallery-plot

   High resolution plots can be downloaded using the links below each plot

