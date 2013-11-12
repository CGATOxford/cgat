=============================
Differential methylation
=============================

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

.. report:: DifferentialMethylation.TrackerDMRSummary
   :render: table

   Table of with results of differential enrichment analysis.

.. ifconfig:: "deseq" in METHODS

   .. toctree::
      DESeq.rst

.. ifconfig:: "deseq" in METHODS

   .. toctree::
      EdgeR.rst





