.. _BindingPatterns:

================
Binding patterns
================

This section provides an overview of where intervals are located with
respect to protein coding genes. The plots below count the number of
genes with at least one binding event in any of the various regions
(cds, utr, etc) that constitute a gene. Note that regions have very
different sizes which distorts these plots somewhat.

.. report:: Binding.BindingSummary
   :render: matrix-plot
   :transform-matrix: normalized-col-first
   :max-rows: 0
   :max-cols: 0
   :palette: Blues
   :restrict: r(peak)
   
   Number of intervals in regions around genes. This plot counts the
   overlap only with the peak.

.. report:: Binding.BindingSummary
   :render: table
   :restrict: r(peak)
   
   Number of intervals in regions around genes. This plot counts the
   overlap only with the peak.


.. report:: Binding.BindingSummary
   :render: matrix-plot
   :transform-matrix: normalized-col-first
   :max-rows: 0
   :max-cols: 0
   :palette: Blues
   :exclude: r(peak)
   
   Number of intervals in regions around genes. This plot counts the
   overlap with the full extent of an interval.

.. report:: Binding.BindingSummary
   :render: table
   :exclude: r(peak)
   
   Number of intervals in regions around genes. This plot counts the
   overlap with the full extent of an interval.

.. report:: Binding.BindingSummary
   :render: table

   Number of intervals in regions around genes.


Binding patterns
=================

The location of intervals with respect to transcripts is summarized by
a pattern. A binding pattern is a string of 0's and 1's delineating the absence
or presence of an interval within a certain region:

+--------------------+--------------------+--------------------+
|Position            |Region              |Description         |
+--------------------+--------------------+--------------------+
|1                   |flank5              |10kb flank upstream |
|                    |                    |of UTR              |
+--------------------+--------------------+--------------------+
|2                   |utr5                |5' UTR              |
+--------------------+--------------------+--------------------+
|3                   |cds                 |coding exon         |
+--------------------+--------------------+--------------------+
|4                   |intron              |intronic regions    |
+--------------------+--------------------+--------------------+
|5                   |utr3                |3' UTR              |
+--------------------+--------------------+--------------------+
|6                   |flank3              |10 kb flank         |
|                    |                    |downstream of UTR   |
+--------------------+--------------------+--------------------+

.. report:: Binding.BindingPatterns
   :render: matrix-plot
   :transform-matrix: normalized-row-total
   :force:
   :tf-labels: pattern
   :max-rows: 0
   :max-cols: 0
   :groupby: slice

   Proportion of transcripts with intervals 
   Note that transcripts without any binding have been removed.

.. report:: Binding.BindingPatterns
   :render: matrix
   :force:
   :tf-labels: pattern
   :max-rows: 0
   :max-cols: 0
   :groupby: slice

   Number of transcripts with intervals. Note that transcripts without
   any binding have been removed.


