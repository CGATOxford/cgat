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

   Percent :term:`model coverage` of transcripts in the reference gene set.

.. report:: Reference.TranscriptCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:

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

   Percent :term:`model coverage` of genes in the reference gene set.

.. report:: Reference.GeneCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:

   Percent :term:`model coverage` of genes in the reference gene set.
   This plot is cumulative.

Mean coverage versus maximum coverage
-------------------------------------
The following plot shows the correlation of mean read depth and
maxmimum read depth. The correlation usually breaks down for long
genes.

.. report:: Reference.MeanVsMaxReadDepth
   :render: scatter-rainbow-plot
   :layout: grid
   :width: 300

   Maxmimum read depth versus mean read depth of :term:`reference` genes.
   Dots are coloured by the log(length) of a :term:`reference` gene.

In contrast, mean and median are usually well correlated:

.. report:: Reference.MeanVsMedianReadDepth
   :render: scatter-rainbow-plot
   :layout: grid
   :width: 300

   Maxmimum read depth versus median read depth of :term:`reference` genes.
   Dots are coloured by the log(length) of a :term:`reference` gene.


