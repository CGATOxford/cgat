========
Coverage
========

Coverage of Target Regions
==========================

The following table present an overview of the sequencing depth acheived for each of the target regions for each :term:`track`.

.. report:: Coverage.CoverageSummary
   :render: table
   :tracks: coverage_stats
   :slices: features, covered, mean_cov, median_cov

   Coverage Summary

.. report:: Coverage.CoveragePlot
   :render: line-plot
   :transform: histogram
   :tf-range: 0,100,1
   :tracks: coverage_stats
   :slices: cov_mean

   Coverage Summary

.. report:: Coverage.CoveragePlot2
   :render: line-plot
   :as-lines:
   :yrange: 0,

   Coverage Summary
