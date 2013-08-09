.. _assembly:


=====================
Meta genome assembly
=====================

This page summarises the reults of the metagenome assembly. In these steps the short-read data are
assembled into contigs for downstream gene estimation and so forth. Common assembly statistics are 
used to charcterise the success of the assembly.



Summary
========

Below is a summary of the assembly results for each sample


.. report:: Assembly.Summary
   :render: table
   

   Summary of contig assembly



Contig lengths
===============

Here we describe the lengths of the contigs. Generally a good assembly will contain contigs > 500bp, which
will enhance the ability to accurately predict open reading frames.



.. report:: Assembly.LengthDistribution
   :render: line-plot
   :transform: histogram
   

   Contig length distribution




