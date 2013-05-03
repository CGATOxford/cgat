=========================
Gene list analysis
=========================

This page summarizes the results of gene list analyses.

.. report:: Genelists.GenelistSummary
   :render: table

   Number of categories called significant in each of the gene list
   analyses

.. report:: Genelists.GenelistSignificantCategories
   :render: matrix-plot
   :transform-matrix: correspondence-analysis
   :force:
   :format: %6.4f
   :palette: RdBu
   :zrange: -2,2
   :layout: column-3
   :width: 300

   Categories enriched/depleted in the various analysis. This set
   is filtered for a fdr. Shown for each category is the log10 fold 
   enrichment. Categories not significant are set to 0.
   

Quality control
===============

Distribution of P-values
-------------------------

The plots below show the distribution of P-Values in various
analyses. For the q-value estimation to work, P-Values should
be almost uniformly distributed with a long straight part in the
middle of the plot.

An excess of small P-Values will indicate that there are not many
true null hypotheses. The q-value-method requires a certain proportion
of true null hypotheses in order to estimate the FDR.

.. report:: Genelists.GenelistPValues
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: arange(0,1.005,0.005)

   Pvalue distribution of tests within each gene list analysis
