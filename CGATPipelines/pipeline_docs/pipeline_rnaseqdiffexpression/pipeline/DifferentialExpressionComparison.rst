=======================
Overlap between methods
=======================

This section compares the results of various methods to predict differential 
expression.

All comparisons are on the gene level, which might not be the best choice
for all methods, for example :term:`cuffdiff`.

Shared genes
============

The following table lists the number of features called differentially
expressed that are common to all methods.

.. report:: DifferentialExpression.DifferentialExpressionOverlap
   :render: table

   Number of consistent predictions between different methods.

.. The following plots all take a while, disabled for now until 
.. there is grouping on a per-design level.

..
   Fold change
   ===========

   .. report:: DifferentialExpression.DifferentialExpressionCorrelationFoldChangeCuffdiffDeseq
      :render: r-smooth-scatter-plot
      :layout: column-3
      :width: 300

      Scatter-plots of fold change estimates Cuffdiff vs. DESeq

   .. report:: DifferentialExpression.DifferentialExpressionCorrelationFoldChangeCuffdiffEdger
      :render: r-smooth-scatter-plot
      :layout: column-3
      :width: 300

      Scatter-plots of fold change estimates Cuffdiff vs. EdgeR

   .. report:: DifferentialExpression.DifferentialExpressionCorrelationFoldChangeDeseqEdger
      :render: r-smooth-scatter-plot
      :layout: column-3
      :width: 300

      Scatter-plots of fold change estimates DESeq vs. EdgeR

   P-values
   ===========

   .. report:: DifferentialExpression.DifferentialExpressionCorrelationPValueCuffdiffDeseq
      :render: r-smooth-scatter-plot
      :layout: column-3

      Scatter-plots of fold change estimates Cuffdiff vs. DESeq


   .. report:: DifferentialExpression.DifferentialExpressionCorrelationPValueCuffdiffEdger
      :render: r-smooth-scatter-plot
      :layout: column-3
      :width: 300

      Scatter-plots of fold change estimates Cuffdiff vs. EdgeR

   .. report:: DifferentialExpression.DifferentialExpressionCorrelationPValueDeseqEdger
      :render: r-smooth-scatter-plot
      :layout: column-3
      :width: 300

      Scatter-plots of fold change estimates DESeq vs. EdgeR





