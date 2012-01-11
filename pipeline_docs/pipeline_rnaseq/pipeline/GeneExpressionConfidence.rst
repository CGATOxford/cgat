.. _FPKM confidence intervals:

=========================
FPKM confidence intervals
=========================
   
This section summarizes the distribution of FPKM confidence
intervals in each of the experiments. Only confidence intervals
for expressed genes, isoforms, etc are shown 
(FPKM > :param:`Expression.ExpressionFPKMConfidence.min_fpkm`).

.. report:: Expression.ExpressionFPKMConfidence
   :render: line-plot
   :as-lines:
   :transform: histogram
   :tf-range: ,,0.1 
   :logscale: x
   :groupby: track
   :layout: column-2
   :width: 300

   Distribution of FPKM values

.. report:: Expression.ExpressionFPKMConfidence
   :render: table
   :transform: stats
   :groupby: track

   Summary of expression FPKM values

.. Correlation between FPKM and confidence intervals
.. =================================================

.. .. report:: Expression.ExpressionFPKMConfidenceCorrelation
..    :render: scatter-plot
..    :logscale: xy
..    :groupby: track
..    :layout: column-2
..    :width: 300

..    FPKM confidence values against FPKM.
