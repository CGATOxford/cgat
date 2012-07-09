==================
Exported data sets
==================

This page describes data sets the pipeline exports per default.

All these files can be found in the project download location.

The download directory contains the following directories:

   genesets
	the various gene sets

   classification
        tables with transcript classes
  
   differential_expression
        differential expression results

Genesets
========

The directory :file:`genesets` contains various gene sets
in :term:`gtf` format. 

abinitio.gtf.gz
   The :term:`abinitio gene set`.

lincrna.gtf.gz
   The :term:`lincrna` gene set.

novel.gtf.gz
   The :term:`novel` gene set.

Tables
======

Tables are in tab-separated format.

Classification of transcripts
-----------------------------

lincrna.class.tsv.gz
	table containing a classification of :term:`lincrna` transcripts

abinitio.class.tsv.gz
	table containing a classification of :term:`abinitio` transcripts

novel.class.tsv.gz
	table containing a classification of :term:`novel` transcripts

Differentially expressed transcripts
------------------------------------

This directory contains the results of the :term:`cuffdiff` tests 
for differential gene expression.







   
