==================
FPKM distributions
==================
   
This section summarizes the distribution of FPKM-values in each of
the experiments. 

Distribution of FPKM values
===========================

.. report:: Expression.ExpressionFPKM
   :render: line-plot
   :as-lines:
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,cumulative
   :groupby: track
   :layout: column-2

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
