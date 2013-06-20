==================
FPKM distributions
==================
   
This section summarizes the distribution of FPKM-values in each of
the experiments. 

Distribution of FPKM values in coding transcripts
=================================================

The following plots show the cumulative distribution of FPKM values.
To facilitate comparison between datasets, the FPKM values are 
normalized by the maximum FPKM per experiment. The more inverse L-shaped
the curve, the more the data set is dominated by highly expressed genes.

.. report:: Expression.ExpressionNormalizedFPKM
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :groupby: track
   :tracks: r(refcoding.*gene)
   :layout: column-3
   :width: 300
   :xrange: ,0.5

   Distribution of normalized :term:`FPKM` values in each experiment.

The following plots show the distribution of unnormalized FPKM values.

.. report:: Expression.ExpressionFPKM
   :render: line-plot
   :as-lines:
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :groupby: track
   :layout: column-3
   :tracks: r(refcoding.*gene)
   :width: 300
   :tf-bins: log-100

   Distribution of FPKM values

Distribution of FPKM values in non-coding sequence
=================================================

The following plots show the cumulative distribution of FPKM values.
To facilitate comparison between datasets, the FPKM values are 
normalized by the maximum FPKM per experiment. The more inverse L-shaped
the curve, the more the data set is dominated by highly expressed genes.

.. report:: Expression.ExpressionNormalizedFPKM
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :as-lines:
   :groupby: track
   :tracks: r(refnoncoding.*gene)
   :layout: column-3
   :width: 300
   :xrange: ,0.5

   Distribution of normalized :term:`FPKM` values in each experiment.

The following plots show the distribution of unnormalized FPKM values.

.. report:: Expression.ExpressionFPKM
   :render: line-plot
   :as-lines:
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :groupby: track
   :tracks: r(refcoding.*gene)
   :layout: column-3
   :width: 300
   :tf-bins: log-100

   Distribution of FPKM values

Summary stats of FPKM values
============================

.. report:: Expression.ExpressionFPKM
   :render: table
   :transform: stats
   :groupby: all

   Summary of expression FPKM values

Highest expressed genes
=======================

.. report:: Expression.ExpressionHighestExpressedGenes
   :render: table
   :force:
   :groupby: track

   The ten highest expressed genes in each data set.
