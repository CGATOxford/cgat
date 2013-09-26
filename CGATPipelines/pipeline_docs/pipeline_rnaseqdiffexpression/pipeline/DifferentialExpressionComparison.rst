=======================
Overlap between methods
=======================

This section compares the results of various methods to predict differential 
expression.

All comparisons are on the gene level, which might not be the best choice
for all methods, for example :term:`cuffdiff`.

Shared genes
============

The following table lists the number of shared predictions in each
gene set.

.. report:: DifferentialExpression.DifferentialExpressionOverlap
   :render: table
   :groupby: all

   Number of shared predictions between different methods.

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





