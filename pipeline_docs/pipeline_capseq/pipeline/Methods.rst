=======
Methods
=======

Data
====

The pipeline detects files with the extension `*.fastq.gz` in the working directory.  
Files are named using the following convention:

Sample-condition-replicate.fastq.gz

ChiP samples have the name of the immunoprecipitated factor as the condition. 
Control samples have 'input' as the condition.

Genome Alignment
================

Reads were aligned to the latest version of the appropriate reference genome using Bowtie (v0.12.7). 
Only uniquely mapping reads, with a maximum of two mismatches across the entire read length were used. 
Reads mapping to exactly the same genomic location (duplicates) were removed using Picard (v1.40). 
Alignment statistics were gathering using both bam2stats.py and Picard.

Peak Finding
============

Binding intervals were identified using MACS (v14) using a bandwith of 300 and an mfold range of 4-30.
Binding intervals were filtered by a q-value (FDR) of 0.01. 

Genomic Annotation
==================

Binding intervals were annotated with respect to the latest Ensembl genome annotation for the appropriate species.


Genomic Association
===================

Binding intervals from different samples, relevant ChIP-seq datasets and CpG island prediction algorithms were 
compared using the Genomic Association Tester (GAT). http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/gat/contents.html.


