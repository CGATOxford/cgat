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
   :glob: counts.dir/counts_stats_heatmap.svg
   :width: 200

   Heatmap of overall sample similarity based on all samples
   Clustered using a correlation distance.

Clustering within designs
-------------------------

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: counts.dir/design*_stats_heatmap.svg
   :width: 200
   :layout: column-2

   Heatmap of overall sample similarity for various
   experimental designs. Clustered using a correlation distance.

Ungrouped heatmaps
------------------

.. report:: MedipReport.TagCountsCorrelations
   :render: matrix-plot
   :groupby: all
   :colorbar-format: %5.2f
   :zrange: 0.9,1.0
   :width: 200
   :layout: column-2

   Ungrouped heatmaps of sample similarity. Shown are heatmaps for
   completed data and various experimental designs.

.. _TagsCounts:

Tag counts in windows
=====================

The following tables report how many tag counts are reported
in total across each window.

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

.. _TagsWindows:

Windows
=======

Window statistics shows the distribution of window sizes used in the analysis.

.. report:: MedipReport.WindowsSummary
   :render: table
   :force:

   Number and size of windows

.. report:: MedipReport.WindowsSizes
   :render: line-plot
   :logscale: xy
   :as-lines:

   Distribution of tile size

.. _TagsDuplicates:

Duplicate statistics
====================

Number and percentage of duplicate pairs removed before tag counting.

.. report:: MedipReport.PicardDuplicatesMetrics
   :render: table
   :force:

   Duplication metrics

.. report:: MedipReport.PicardDuplicatesHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,

   Histogram of duplication levels

