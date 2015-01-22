======================================
Correlations, heatmaps and clusterings
======================================
The plots on this page show global comparisons and contrasts between
the samples in the data set. Correlations are shown between pairs of
treatments per tissue. Below this are heatmaps and hierachical
clusterings from spearman's rank correlations between all pairs of
samples. All results are split into CpGs inside and outside CpG islands.


Global correlation
==================
The following plots are correlations between treatment groups for all
CpGs, split by CpG island / Non-CpG island.

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: plots.dir/*global_correlation*.png	    

   High resolution plots can be downloaded using the links below each plot


Spearman's correlation heatmaps
===============================

The following plots are heatmaps of spearman's rho values for all
possible pairwise comparisons. The samples are ordered by treatment
group and are not clustered. For clustering see the plots in the next
section. Heatmaps are split by CpG island / Non-CpG island.

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: plots.dir/*spearmans_heatmap_*.png
 
   High resolution plots can be downloaded using the links below each plot


Hierachical clustering
======================

The following plots are hierarchical clusterings of the samples based
on pairwise spearman's rank correlations. Clusterings are split by are
split by CpG island / Non-CpG island. Distance between samples is
represented by height on the y-axis. The clustering has been
bootstrapped to estimate how strongly the data support the
cluster. The numbers in green are the bootstrap probabilities.

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: plots.dir/*hierarchical_clustering*.png 

   High resolution plots can be downloaded using the links below each plot

