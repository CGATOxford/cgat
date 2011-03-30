==================
Reference coverage
==================

Coverage per gene
=================

The following plots examine the coverage of reference gene model 
by transcript models.

.. report:: Reference.CoverageByTranscripts
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: arange(0,101,1.0)
   :as-lines:

   Percent overlap of genes in the ENSEMBL gene set by transcript models.

.. report:: Reference.OverlapByTranscripts
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: arange(0,101,1.0)
   :as-lines:

   Percent overlap of genes in the ENSEMBL gene set by transcript models.
   This plot is cumulative.

.. _TableReferenceCoverage:

.. report:: Reference.CoverageStats
   :render: table

   Statistics of the overlap of reference genes by transcript models.

Coverage versus gene length
===========================

The following plots correlate three measures the relative coverage of an ENSEMBL gene model (1-100%)
with read depth or statistical coverage of an ENSEMBL gene model. The numbers do not take into account 
of multiple transcripts mapping to an ENSEMBL gene, but only take the best match into account.

.. report:: Reference.LengthVsReadDepthByLength
   :render: scatter-rainbow-plot

   Relative coverage of known gene (1-100%) versus the read depth.
   Dots are colored by the length of the known gene (in log scale).

.. report:: Reference.CoverageVsLengthByReadDepth
   :render: scatter-rainbow-plot

   Plot the absolute coverage of known gene set versus its length.
   Dots are colored by read depth.

Mean coverage versus maximum coverage
=====================================

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

