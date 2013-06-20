======================
Splice site validation
======================

.. report:: BenchmarkReport.ExonValidationSummary
   :render: interleaved-bar-plot
   :slices: spliced

   Number of spliced reads

Spliced reads
-------------

.. report:: BenchmarkReport.ExonValidationSummary
   :render: table

   Proportion of spliced reads overlapping reference exons

.. report:: BenchmarkReport.ExonValidationSummary
   :render: interleaved-bar-plot
   :slices: spliced_nooverlap,spliced_halfoverlap,spliced_bothoverlap

   Proportion of spliced reads overlapping reference exons

.. report:: BenchmarkReport.ExonValidationSummary
   :render: pie-plot
   :slices: spliced_nooverlap,spliced_halfoverlap,spliced_bothoverlap
   :width: 300
   :layout: column-3

   Proportion of spliced reads overlapping reference exons

.. report:: BenchmarkReport.ExonValidationSummary
   :render: matrix
   :slices: spliced_nooverlap,spliced_halfoverlap,spliced_bothoverlap
   :transform-matrix: normalized-row-total

   Proportion of spliced reads overlapping reference exons

.. report:: BenchmarkReport.ExonValidationSummary
   :render: pie-plot
   :slices: spliced_inexact,spliced_exact
   :width: 300
   :layout: column-3

   Proportion of splice sites that are exact.

.. report:: BenchmarkReport.ExonValidationSummary
   :render: matrix
   :slices: spliced_inexact,spliced_exact
   :transform-matrix: normalized-row-total

   Proportion of splice sites that are exact.

Overrun in unspliced reads
--------------------------

.. report:: BenchmarkReport.ExonValidationSummary
   :render: interleaved-bar-plot
   :slices: unspliced_overlap,unspliced_overrun

   Overrun in unspliced reads

.. report:: BenchmarkReport.ExonValidationSummary
   :render: interleaved-bar-plot
   :transform-matrix: normalized-row-max
   :slices: unspliced_overlap,unspliced_overrun

   Overrun in unspliced reads

.. report:: BenchmarkReport.ExonValidationSummary
   :render: matrix
   :transform-matrix: normalized-row-max
   :slices: unspliced_overlap,unspliced_overrun

   Overrun in unspliced reads

