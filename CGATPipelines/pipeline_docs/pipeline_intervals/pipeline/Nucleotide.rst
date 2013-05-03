=============
Nucleotide
=============

The GC composition of all intervals is examined

Histograms
==========

.. report:: Nucleotide.pGC
   :render: line-plot
   :transform: histogram
   :tf-bins: arange(0.2,0.8,0.03)
   :as-lines:


   Histogram of interval GC contents


GC vs average coverage
======================

.. report:: Nucleotide.GCvsCounts
    :render: scatter-plot
    :groupby: track

    Scatter plots of interval GC vs average no.counts


Gene profiles of genes with and without interval overlap
========================================================

.. report:: GeneProfile.MultiProfilePlot
   :render: line-plot
   :groupby: track
   :as-lines: 

     Read coverage of gene models either overlapping or not-overlapping intervals

GC content of 100bp windows surrounding the TSS of genes
========================================================

.. report:: Nucleotide.TSS
   :render: box-plot
   :groupby: tracks
   :mpl-rc: figure.subplot.bottom=0.2

   Boxplots of the GC content of 100bp windows of gene TSS. Genes are split into those overlap the intervals and those that do not (any of the gene model can intersect with an interval to qualify it as overlapping).


.. report:: Nucleotide.TSS
   :render: line-plot
   :transform: histogram
   :tf-bins: arange(0.2,0.8,0.03)
   :as-lines:

   Histogram of interval GC contents


