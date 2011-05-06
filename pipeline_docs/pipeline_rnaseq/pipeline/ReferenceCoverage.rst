===================
Reference gene set
===================

This section looks at the reference gene set and the reads.

Coverage
========

Gene coverage
--------------

The following plots examine the :term:`model coverage` of transcripts and genes
in the reference gene set by reads. Note that only coverage of sense reads is
reported.

.. report:: Reference.TranscriptCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:

   Percent :term:`model coverage` of genes in the reference gene set.

.. report:: Reference.GeneCoverage
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: arange(0,101,1.0)
   :as-lines:

   Percent :term:`model coverage` of genes in the reference gene set.
   This plot is cumulative.

Coverage versus gene length
---------------------------

The following plots correlate three measures the relative coverage of a reference gene model (1-100%)
with read depth or :term:`depth coverage` of the reference gene model. 

.. report:: Reference.CoverageVsLengthByReadDepth
   :render: scatter-rainbow-plot

   Plot the absolute coverage of known gene set versus its length.
   Dots are colored by read depth.

Mean coverage versus maximum coverage
-------------------------------------
The following plot shows the correlation of mean read depth and
maxmimum read depth. The correlation usually breaks down for long
genes.

.. report:: Reference.MeanVsMaxReadDepth
   :render: scatter-rainbow-plot
   :layout: grid

   Maxmimum read depth versus mean read depth of :term:`reference` genes. 
   Dots are coloured by the log(length) of a :term:`reference` gene.

In contrast, mean and median are usually well correlated:

.. report:: Reference.MeanVsMedianReadDepth
   :render: scatter-rainbow-plot
   :layout: grid

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

   Directionality of reads within transcript models.

