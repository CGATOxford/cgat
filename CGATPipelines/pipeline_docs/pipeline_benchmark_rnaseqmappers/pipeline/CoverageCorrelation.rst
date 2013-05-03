==========================
Correlation of read depths
==========================

Exons
=====

.. report:: BenchmarkReport.CoverageCorrelationExons
   :render: matrix-plot
   :transform: pearson,select
   :tf-fields: coefficient
   :transform-matrix: correspondence-analysis
   :force:
   :zrange: 0.8,1
   :colorbar-format: %6.4f

   Correlation coefficients: number of reads aligned
   to genes between different methods

.. report:: BenchmarkReport.CoverageCorrelationExons
   :render: matrix
   :transform: pearson,select
   :logscale: xy
   :tf-fields: coefficient
   :force:
   :format: %6.4f

   Correlation coefficients: number of reads aligned
   to genes between different methods
