=======================
Hypergeometric analysis
=======================

This section lists the results from hypergeometric enrichment analysis
on various gene lists. See :ref:`genelists`.

Hypergeometric pathway analysis is the most basic analysis of pathway
enrichment. The hypergeometric test evaluates if the number of genes
in a gene list that occur in a given pathway are more or less frequent
than expected by chance.

The tests do not account for any knowledge about the interaction
between genes in a pathway (activator, repressor). Neither do the tests 
correct for any possible biases that occur through the derivation of
the genelists or pathways.

.. report:: Hypergeometric.HypergeometricResults
   :render: table
   :groupby: none
   :force:
                                             
   Tabular results of hypergeometric enrichment analysis.

..
   .. report:: Hypergeometric.FoldMatrix
      :render: r-heatmap-plot
      :extra-formats: svg

      GO enrichment of various gene sets.

..
   Clustered analysis
   ==================

   GO results have been filtered removing semantically similar GO terms.
   The size of a dot indicates the significance of a result.

   .. report:: Tracker.TrackerImages
      :render: gallery-plot
      :glob: hypergeometric.dir/*/*.svg

      Semantic clustering of GO results

..
   GO analyses
   ===========

   The following tables list the top 20 significant results for various
   gene list definitions.

   Promotor definitions
   --------------------

   Definitions of genelists by promotor occupancy. The promotor is a
   region of -2kb/+0.5kb around the TSS.

   This definition is valid for each FDR separately.

   .. report:: Genelists.GOAnalysisComparisonPromotor
      :render: table

      Using Promotor definitions

   .. report:: Genelists.GOAnalysisComparisonPromotor
      :render: interleaved-bar-plot
      :transform: tolabels
      :tracks: biol_process
      :xrange: 0,5
      :switch:
      :orientation: horizontal
      :tf-labels: description

      Using Promotor definitions

