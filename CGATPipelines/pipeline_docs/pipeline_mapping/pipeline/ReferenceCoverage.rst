.. _Reference Coverage:

===================
Reference coverage
===================

This section looks at the reference gene set and the reads.

Coverage
========

Transcript coverage
--------------------

The following plots examine the :term:`model coverage` of transcripts
in the reference gene set by reads. Coverage is computed for each
transcript separately. Reads will contribute coverage to multiple
transcripts if the read overlaps several transcripts.

Only sense reads are reported. For protocols that do not
preserve strand only half the reads will be used to calculate coverage.

.. report:: Reference.TranscriptCoverage
   :render: table
   :transform: stats
   :force:

   Percent :term:`model coverage` of transcripts in the reference gene set.

.. report:: Reference.TranscriptCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:
   :layout: column-3

   Percent :term:`model coverage` of transcripts in the reference gene set.

Gene coverage
--------------------

The following plots examine the :term:`model coverage` of genes
in the reference gene set by reads. Coverage is computed for each
transcript separately. Reads will contribute coverage to multiple
transcripts if the read overlaps several transcripts. The coverage 
per gene is the maximum coverage of a transcript.

Only sense reads are reported. For protocols that do not
preserve strand only half the reads will be used to calculate coverage.

.. report:: Reference.GeneCoverage
   :render: table
   :transform: stats
   :force:

   Percent :term:`model coverage` of genes in the reference gene set.

.. report:: Reference.GeneCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:
   :layout: column-3

   Percent :term:`model coverage` of genes in the reference gene set.
   This plot is cumulative.

Directionality
==============

This section looks at the directionality of reads within transcript models.
Libraries without strand information should have a peak at about 1.0.

.. report:: Reference.ReadDirectionality
   :render: line-plot
   :transform: histogram
   :slices: transcript
   :logscale: x
   :tf-aggregate: normalized-total
   :tf-range: ,,0.1
   :groupby: slice
   :as-lines:
   :layout: column-3
   :width: 300

   Directionality of reads within transcript models.

