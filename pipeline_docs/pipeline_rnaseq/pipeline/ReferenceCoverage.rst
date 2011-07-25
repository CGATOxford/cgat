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

Full length transcript models
-----------------------------

The following plot examines the proportion of full length transcript
models that can theoretically be built at a certain expression level.
Each plot contains several lines corresponding to what constitutes a
full length transcript (>50% of transcript is covered with reads up to
>90% of transcript is covered with reads).

.. report:: Reference.CoverageVsFPKM
   :render: line-plot
   :logscale: x
   :groupby: track
   :as-lines:
   :width: 200
   :layout: column-3

   Proportion of full-length transcripts versus FPKM. 
   
These curves converge to the number of full length transcripts, which
is different depending on what percentage of a transcript is considered to be full-length.

Coverage versus gene length
---------------------------

The following plots correlate three measures the relative coverage of a reference gene model (1-100%)
with the expression level of a transcript.

.. report:: Reference.CoverageVsLengthByReadDepth
   :render: scatter-rainbow-plot
   :layout: column-3
   :width: 300
   :mpl-rc: lines.markersize=2

   Plot the coverage of a transcript versus its expression level.
   Dots are colored by transcript length.

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

Directionality
==============

This section looks at the directionality of reads within transcript models.
Libraries without strand information should have a peak at about 1.0.

.. report:: Reference.ReadDirectionality
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :tf-range: ,,0.1
   :groupby: slice
   :as-lines:
   :layout: column-3
   :width: 300

   Directionality of reads within transcript models.

