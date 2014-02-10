===================
Cuffdiff Results
===================

Differential expression
=============================

Summary
-------

The following table summarizes the results of the test for differential expression.
It lists for every combination of object (gene, isoform, ...), geneset and pair of tracks
the number of tests that passed, failed, etc.

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

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: table

   Table of Cuffdiff differential expression results

.. report:: DifferentialExpression.TrackerDESummaryCuffdiff
   :render: interleaved-bar-plot
   :force:
   :slices: significant

   Number of significant tests

Diagnostic plots
----------------

.. report:: Tracker.TrackerImages 
   :render: gallery-plot    
   :glob: cuffdiff.dir/*.png
   
   Diagnostic plots from the Cuffdiff analysis
