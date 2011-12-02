=============================
Differential methylation
=============================

The following sections look at the results of predicting
regions of differential methylation.

The following table summarizes the results of the test for
differential methylation using :term:`DESeq` and :term:`EdgeR`.

+----------------------------------------+----------------------------------------+
|*Column*                                |*Content*                               |
+----------------------------------------+----------------------------------------+
|significant                             |Number of significant tests             |
+----------------------------------------+----------------------------------------+
|tested                                  |Number of tests in total                |
+----------------------------------------+----------------------------------------+
|status_OK                               |Number of successful tests              |
+----------------------------------------+----------------------------------------+
|twofold                                 |Number of significant tests, where the  |
|                                        |expression level differs by at least a  |
|                                        |factor of 2.                            |
+----------------------------------------+----------------------------------------+

.. report:: DifferentialMethylation.TrackerDMRSummaryDESeq
   :render: table

   Table of DESeq differential expression results

.. report:: DifferentialMethylation.TrackerDMRSummaryEdgeR
   :render: table

   Table of EdgeR differential expression results

.. toctree::
   :maxdepth: 2

   Tiling.rst
      

