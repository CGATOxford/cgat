==========================
Transcript and gene models
==========================

This section contains summaries and statistics about the transcript and gene-model
build process within the RNASeq pipeline.

Gene set
========

The following sections summarize the contents from the final 
gene sets. 

.. toctree::
   :maxdepth: 2

   GeneModelsFinalSummary.rst
   GeneModelsFinalMappability.rst
   ReferenceCoverage.rst

Intermediate results
====================

The sections below track intermediate results in the various stages of the 
gene model building process.

.. toctree::
   :maxdepth: 2

   GeneModelsBenchmark.rst
   GeneModelsTranscripts.rst
   GeneModelsReproducibility.rst
   GeneModelsSummary.rst

.. _CuffCompare:

Cuffcompare overview
====================

:term:`cuffcompare` compares the predicted transcripts from one or more experiments amongst each other. 
Similar transcripts are aggregated as :term:`transfrags` and designated a shared identifier permitting
similar transcripts to be tracked across experiments. Furthermore, each transfrag is compared to
the reference gene set and classified into the following categories:

+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|Priority|Code|Description                                                                                                                                              |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|1       | =  | Complete match of intron chain                                                                                                                          |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|2       | c  | Contained                                                                                                                                               |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|3       | j  | Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript                                                |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|4       | e  | Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron indicating a possible pre-mRNA fragment.                    |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|5       | i  | A transfrag falling entirely within a reference intron                                                                                                  |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|6       | o  | Generic exonic overlap with a reference transcript                                                                                                      |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|7       | p  | Possible polymerase run-on fragment (within 2Kbases of a reference transcript)                                                                          |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|8       | r  | Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|9       | u  | Unknown intergenic transcript                                                                                                                           |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|10      | x  | Exonic overlap with reference on the opposite strand                                                                                                    |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|11      | s  | An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)                                       |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|12      | .  | indicates multiple classifications                                                                                                                      |
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+
|        | *  | Any of the above - used to display summaries on all transfrags                                                                                          |
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+

