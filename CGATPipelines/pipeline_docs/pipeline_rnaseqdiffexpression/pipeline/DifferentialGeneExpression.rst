============================
Differential gene expression
============================

The following sections reports the number of windows that contain a
significantly different number of reads.


+--------------------+------------------------------------------------------------+
|*Column*            |*Content*                                                   |
+--------------------+------------------------------------------------------------+
|method              |Method used for differential analysis                       |
+--------------------+------------------------------------------------------------+
|treatment           |Treatment name                                              |
+--------------------+------------------------------------------------------------+
|control             |Control name                                                |
+--------------------+------------------------------------------------------------+
|up                  |Total number of up-regulated windows                        |
+--------------------+------------------------------------------------------------+
|down                |Total number of down-regulated windows                      |
+--------------------+------------------------------------------------------------+
|tested              |Number of tested windows in total                           |
+--------------------+------------------------------------------------------------+
|significant         |Number of significant windows                               |
+--------------------+------------------------------------------------------------+
|fold2               |Number of significant windows with > 2-fold difference      |
+--------------------+------------------------------------------------------------+
|significant_down    |Number of down-regulated windows called significant         |
+--------------------+------------------------------------------------------------+
|significant_up      |Number of up-regulated windows called significant           |
+--------------------+------------------------------------------------------------+
|l2fold_down         |Number of down-regulated windows called significant with >  |
|                    |2-fold difference                                           |
+--------------------+------------------------------------------------------------+
|l2fold_up           |Number of up-regulated windows called significant with >    |
|                    |2-fold difference                                           |
+--------------------+------------------------------------------------------------+

.. report:: DifferentialExpression.DESummary
   :render: table

   Table with results of differential enrichment analysis.

More detailed results for each method are below.

.. ifconfig:: "deseq" in PARAMS['methods']

   .. toctree::
      DESeqResults.rst

.. ifconfig:: "deseq" in PARAMS['methods']

   .. toctree::
      EdgeRResults.rst

.. ifconfig:: "cuffdiff" in PARAMS['methods']

   .. toctree::
      CuffdiffResults.rst

   .. toctree::
      DifferentialExpressionComparison.rst

