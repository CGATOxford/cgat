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
|geneset             |Geneset being analysed                                      |
+--------------------+------------------------------------------------------------+
|level               |Object being assessed (transcript/gene)                     |
+--------------------+------------------------------------------------------------+
|tested              |Number of tests run                                         |
+--------------------+------------------------------------------------------------+
|significant         |Number of significant tests                                 |
+--------------------+------------------------------------------------------------+
|twofold             |Number of tests with a two-fold change difference           |
+--------------------+------------------------------------------------------------+

.. report:: DifferentialExpression.DESummary
   :render: table
   :slices: significant,twofold,tested

   Table with results of differential enrichment analysis.

.. Turned off - needs proper grouping

   .. report:: DifferentialExpression.DESummary
      :render: interleaved-bar-plot
      :slices: significant

      Number of significant number of DE genes

   .. report:: DifferentialExpression.DESummary
      :render: interleaved-bar-plot
      :slices: twofold

      Number of twofold number of DE genes

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

