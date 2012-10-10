


=========
Summary
=========

This section of the report summarises the prediction of lncRNA from ab initio built transcripts.
The summary reports the results of the filtering from ab initio assembled lncRNA;


Step 1 - abinitio lncRNA set
==============================

lncRNA are predicted from an ab initio assembly of transcripts. Transcripts are filtered out that 
overlap ENSEMBL transcripts with the following biotypes:

protein_coding\nprocessed_pseudogene\nunprocessed_pseudogene\nnonsense_mediated_decay\nretained_intron\n

+------------------+------------------+
| Input            | Output           |
|                  |                  |
+------------------+------------------+
| <name>.gtf.gz    | lncrna.gtf.gz    |
|                  |                  |
+------------------+------------------+

They are also filtered to be of length > 200bp

Step 2 - filter abinitio lncRNA set
====================================

During an initial assembly using cufflinks, there are a number of fragmentory and single exon
transcripts. Single exon transcripts (i.e. not supported by splice junctions) are generally considered to
be low confidence lncRNA. Therefore, at this point we filter out any single exon lncRNA. However, we include 
a parameter to specify a set of 'known' lncRNA, as well as a comparison of a non-coding RNA set from ensembl. If our
single exon lncRNA overlap with a single exon lncRNA from the reference sets then the lncRNA is kept (from the
reference). This allows us to retain single exon transcripts that have multiple lines of evidence supporting
their existence.


+--------------------+------------------------+
| Input              | Output                 |
|                    |                        |
+--------------------+------------------------+
| lncrna.gtf.gz      | lncrna_filtered.gtf.gz |
| <reference>.gtf.gz |                        |
+--------------------+------------------------+


Step 3 - coding potential calculation
=====================================

Predicted lncRNA from the above steps should be tested for their poteential to code for proteins. We therefore
run the coding potential calculator on the transcripts in this set (`PMC1933232`_).
If a gene contains a transcript that has a coding potetential then all transcripts for this gene are removed at this
stage.

.. _PMC1933232: http://www.ncbi.nlm.nih.gov/pubmed/17631615


+-------------------------+------------------------+
| Input                   | Output                 |
|                         |                        |
+-------------------------+------------------------+
| lncrna_filtered.gtf.gz  | lncrna_final.gtf.gz    |
|                         |                        |
+-------------------------+------------------------+


Step 4 - classification of lncRNA
==================================

Predicted lncRNA are classified based on their proximity to protein coding loci.

+-------------------------+---------------------------+
| Input                   | Output                    |
|                         |                           |
+-------------------------+---------------------------+
| lncrna_final.gtf.gz     | lncrna_final.class.gtf.gz |
|                         |                           |
+-------------------------+---------------------------+


Overview
========

Below is an overview of the number of lncRNA that remain at each stage of the filtering process.


.. report:: Summary.Summary
   :render: table
   
   Summary of lincRNA identified at each stage


.. report:: Summary.Summary
   :render: interleaved-bar-plot
   
   Summary of lincRNA identified at each stage







