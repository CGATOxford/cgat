=============
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

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: *_counts.dir/*_counts_stats_heatmap.svg
   :width: 200

   Heatmap of overall sample similarity based on all samples
   Clustered using a correlation distance.

Clustering within designs
-------------------------

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: *_counts.dir/design*_stats_heatmap.svg
   :width: 200
   :layout: column-2

   Heatmap of overall sample similarity for various
   experimental designs. Clustered using a correlation distance.

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
   completed data and various experimental designs.

.. _TagsCounts:

Tag counts in features
======================

The following tables report how many tag counts are reported
in total across each feature

.. report:: TagCounting.TagCountsSummaryAll
   :render: table
   :force:

   Number of features with a certain number of 
   tag counts based on all samples.

.. report:: TagCounting.TagCountsSummaryPerDesign
   :render: table
   :force:

   Number of features with a certain number of 
   tag counts for various experimental designs.
