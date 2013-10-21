.. _functionalProfiling:


====================
Functional profiles
====================

Metagenomic data analysis often seeks to answer questions regarding the functional
capabilities of the community. This can be achieved in two ways. Firstly, and most
conveniently (though not neccessarily most accurately) this can be done through alignment
of raw reads to functional annotations. The pipeline attempts to carry out this analysis
by using RPS-BLAST to align reads to a database of conserved protein domains (NCBI CDD). These
are then annotated with their associated Cluster of Orthologous Groups (COG) identifier and
functional categories represented by the reads are counted.


Reads aligning to conserved protein domains
============================================

The number of reads aligning to functional domains are provided. This gives us an idea of
the quality of the sequencing as the majority of reads from bacterial species shouls align 
to a known protein domain.


.. report:: FunctionalProfile.AlignmentCounts
   :render: table


   Summary of the number of reads aligning to protein domains


Functional groups repersented by reads
=======================================

It is of interest to know which functional groups are most heavily represented in the
community. Below is a summary of functional groups (COGs) that are represented by the
reads.

.. report:: FunctionalProfile.CogCounts
   :render: table

   Summary of functional categories represented by reads


.. report:: FunctionalProfile.CogCounts
   :render: interleaved-bar-plot
   :mpl-rc: figure.figsize=(20,10);legend.fontsize=10
   :legend-location: lower-right

   Summary of functional categories represented by reads

