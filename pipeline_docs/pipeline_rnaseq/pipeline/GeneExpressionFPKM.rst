==================
FPKM distributions
==================
   
This section summarizes the distribution of FPKM-values in each of
the experiments. 

Distribution of FPKM values
===========================

The following plots show the cumulative distribution of FPKM values.
To facilitate comparison between datasets, the FPKM values are 
normalized by the maximum FPKM per experiment. The more inverse L-shaped
the curve, the more the data set is dominated by highly expressed genes.

.. report:: Expression.ExpressionNormalizedFPKM
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :xrange: 0,0.4
   :as-lines:
   :groupby: track
   :tracks: r(refcoding.*gene)

   Distribution of normalized :term:`FPKM` values in each experiment.

The following plots show the distribution of unnormalized FPKM values.

.. report:: Expression.ExpressionFPKM
   :render: line-plot
   :as-lines:
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :groupby: track
   :layout: column-2
   :tracks: r(refcoding.*gene)

   Distribution of FPKM values

Summary stats of FPKM values
============================

.. report:: Expression.ExpressionFPKM
   :render: table
   :transform: stats
   :groupby: track

   Summary of expression FPKM values

Highest expressed genes
=======================

.. report:: Expression.ExpressionHighestExpressedGenes
   :render: table
   :force:
   :groupby: track

   The ten highest expressed genes in each data set.
