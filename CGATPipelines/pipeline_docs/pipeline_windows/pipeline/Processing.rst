=============
Tag counting
=============

Sample clustering
=================

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :tracker: counts.dir/counts_stats_heatmap.svg

   Heatmap of overall sample similarity based on all samples
   Clustered using a correlation distance.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :tracker: counts.dir/design*_stats_heatmap.svg

   Heatmap of overall sample similarity for various
   experimental designs. Clustered using a correlation distance.

.. report:: MedipReport.TagCountsCorrelations
   :render: matrix-plot
   :groupby: all
   :colorbar-format: %5.2f
   :zrange: 0.9,1.0

   Ungrouped heatmaps of sample similarity. Shown are heatmaps for
   completed data and various experimental designs.

Tag counts in windows
=====================

.. report:: MedipReport.TagCountsSummary
   :render: table
   :force:

   Number of windows with a certain number of 
   tag counts based on all samples.

.. report:: MedipReport.TagCountsDesignSummary
   :render: table
   :force:

   Number of windows with a certain number of 
   tag counts for various experimental designs.

Windows
=======

Window statistics

.. report:: MedipReport.WindowsSummary
   :render: table
   :force:

   Number and size of windows

.. report:: MedipReport.WindowsSizes
   :render: line-plot
   :logscale: xy
   :as-lines:

   Distribution of tile size

Duplicate statistics
====================

.. report:: MedipReport.PicardDuplicatesMetrics
   :render: table
   :force:

   Duplication metrics

.. report:: MedipReport.PicardDuplicatesHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,

   Histogram of duplication levels

