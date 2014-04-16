C=============
Tag counting
=============

.. _TagsSampleClustering:

Sample clustering
=================

The following plots show how the various samples cluster according to
the different experimental designs that have been submitted to the
pipeline.

The clustering is performed using a read-count correlation distance.

Overall clustering
------------------

The plots below show the clustering of all samples using a correlation
distance. There is one plot for each combination of gene set and
counting method.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: *_counts.dir/*_counts_stats_heatmap.svg
   :width: 200

   Heatmap of overall sample similarity based on all samples
   Clustered using a correlation distance.

Clustering within designs
-------------------------

The plots below show the clustering of all samples that are part of an
experimental design. There is one plot for each combination of design,
gene set and for each counting method.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: designs.dir/design*_stats_heatmap.svg
   :width: 200
   :layout: column-2

   Heatmap of overall sample similarity for various
   experimental designs. Clustered using a correlation distance.

..
   Ungrouped heatmaps
   ------------------

   .. report:: TagCounting.TagCountsCorrelations
      :render: matrix-plot
      :groupby: all
      :colorbar-format: %5.2f
      :zrange: 0.9,1.0
      :width: 200
      :tight-layout:
      :layout: column-2

      Ungrouped heatmaps of sample similarity. Shown are heatmaps for
      complete data and various experimental designs.

.. _TagsCounts:

Tag counts in features
======================

The following tables report how many tag counts are reported
in total across each feature within the various experimental designs.
For example, ``max_per_row <= 2`` means that for this number of
features at most 2 tags have been observed in any of the samples.

The first table below shows the counts across all samples.

.. report:: TagCounting.TagCountsSummaryAll
   :render: table

   Number of features with a certain number of 
   tag counts based on all samples.

The table below contains one section for each combination of design,
gene set and for each counting method.

.. report:: TagCounting.TagCountsSummaryPerDesign
   :render: table

   Number of features with a certain number of 
   tag counts for various experimental designs.


Feature counts
==============

featureCounts_ is a tool to compate the number of reads or pairs
within a geneset. The following plot shows the feature counts
summary statistics for all tracks:

.. report:: TagCounting.FeatureCountsSummary
   :render: table

   Feature counts summary.

.. report:: TagCounting.FeatureCountsSummary
   :render: stacked-bar-plot
 
   Feature counts summary

   



