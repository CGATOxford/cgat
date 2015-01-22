=======================================
Differentially Methylated Regions - M3D
=======================================

The plots below present the results from the M3D tools developed by
Guido Sanguinetti's group in Edinburgh to identify differentially
methylated regions from changes in methylation profiles at CpG clusters

Please see http://arxiv.org/abs/1410.6677 for details on M3D

In brief, M3D generates a test statistic in the range 0-1 which
represents how similar two methylation profiles are over a CpG cluster
(0 being identical). All possible pairwise comparison are made for
each CpG cluster. Comparisons between replicates are used to generate
a null distribution of the expected M3D test statistics from
biological variability alone. The M3D scores from the inter-group
comparisons are then compared to the null ditribution to generate p-values
empirically. p-values are adjusted to account for multiple testing
using the Benjamini Hochberg adjustment. Adjusted p-values below 0.05
are identified as significantly differentially methylated.

Clusters were defined as a minimum of 10 CpGs with a maximum distance
between neighbouring CpGs of 100 bp.


M3D statistic distibutions
==========================
The following plot shows the distribution of M3D test statistics for
the within replicates and between treatment group comparisons

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: M3D_plots.dir/M3D_*stat_hist*png

   High resolution plots can be downloaded using the links below each plot


M3D statistic distibutions
==========================
The following table shows the number of signficantly differentially
methylated regions detected by M3D out of the total number of clusters.

.. report:: M3D.M3DSummaryTable
   :render: table


M3D power analysis
==================
To assess the power which we have to detect differentially methylated
regions, in-silico clusters were generateed by swapping CpGs between
clusters to generate artificial methylation differences.

For each plot the :term:`size` axis represents the number of
sequential CpGs which were swapped and the :term:`change` axis
represents the magnitude of the difference between the mean
methylation values of the original CpGs and the CpGs which were
swapped in.

For each data point shown, 100 artifical clusters were generated and
the number of signficant DMRs identified was used to calculate power.

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: power.dir/M3D_change_vs_size_power_*.png

   High resolution plots can be downloaded using the links below each plot


.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: power.dir/M3D_powered.png

   The plot above show the mimimum required methylation changes in order to
   reach 0.8 power.

   High resolution plots can be downloaded using the links below each plot


