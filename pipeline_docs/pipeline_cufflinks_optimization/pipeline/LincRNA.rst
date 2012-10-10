
=======================================
Building lincRNA transcripts abinitio
=======================================

Pipeline cufflinks optimization is used to assess the parameter values in cufflinks that are optimal 
for a givekn question and data set. In this instance we are interested in how well we are able to 
build lincRNA transcripts ab initio i.e. we are not relying on a reference transcriptome. In this sense
we are able to build novel transcripts and those that have not been seen before. We summarise the results
of the lincRNA transcript building by assessing:

* The number of multi-exonic lincRNA we build ab initio

* The average number of exons that the multi-exonic lincRNA contain - the higher the better as this is suggestive
of more full length transcripts

* The length of multi exonic transcripts - the longer the better

* The proportion of the total number of lincRNA that are single exons


Number of splice reads contributing to transcripts
===================================================

Provides a rough idea of split-read usage by Cufflinks

.. report:: ReadsContributingToTranscripts.ReadsContributingToTranscriptsSummary
   :render: table
   :slices: percent alignments

   Proportion of reads contributing to transcripts


.. report:: ReadsContributingToTranscripts.ReadsContributingToTranscriptsSummary
   :render: interleaved-bar-plot
   :slices: percent alignments

   Proportion of reads contributing to transcripts


Number of multi-exonic and single exon lincRNA
================================================

.. report:: LincRNAStats.SingleAndMultiExonCounts
   :render: table
   

   Number of single and multi-exon lincRNA



.. report:: LincRNAStats.SingleAndMultiExonCounts
   :render: interleaved-bar-plot
   

   Number of single and multi-exon lincRNA


Number of exons in multi exonic lincRNA
=========================================

.. report:: LincRNAStats.LincRNAExonSummary
   :render: table
   

   Number of exons in multi-exon lincRNA


.. report:: LincRNAStats.LincRNAExonSummary
   :render: interleaved-bar-plot
   

   Number of exons in multi-exon lincRNA



Length of multi-exonic lincRNA
===============================


.. report:: LincRNAStats.LincRNALengthSummary
   :render: table
   

   Average length of multi-exon lincRNA transcripts 


.. report:: LincRNAStats.LincRNALengthSummary
   :render: interleaved-bar-plot
   

   Average length of multi-exon lincRNA transcripts 


Proportion of total lincRNA that are single exons
==================================================


.. report:: LincRNAStats.ProportionSingleExon
   :render: table
   

   Proportion of total lincRNA that are single exon 



.. report:: LincRNAStats.ProportionSingleExon
   :render: interleaved-bar-plot
   

   Proportion of total lincRNA that are single exon 
