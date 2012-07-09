.. _Mapping:

=======
Mapping
=======

This section summarizes results of the read mapping process.

Overview
========

The following table summarizes the number of reads input, mapped and
filtered. 

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|reads_in                               |Number of reads submitted                         |
+---------------------------------------+--------------------------------------------------+
|reads_out                              |Number of reads removed by tophat                 |
+---------------------------------------+--------------------------------------------------+
|reads_mapped                           |Number of reads mapped (after filtering)          |
+---------------------------------------+--------------------------------------------------+
|reads_norna                            |Number of reads not mapping to repetitive RNA     |
+---------------------------------------+--------------------------------------------------+
|reads_norna_unique_alignments          |Number of reads with unique alignments            |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_in,reads_out,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_in,reads_out,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping Results

Detailed results
================

Below are detailed results from the mapping process. It reports summary statistics from various
:term:`bam` formatted files that the pipeline creates.

.. toctree::
   :maxdepth: 2

   MappingProcess.rst
   MappingContext.rst
   MappingAlignmentStatistics.rst

There are various alignment files created in the mapping process:

+---------------------------------------+--------------------------------------------------+
|*Filename*                             |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.accepted.bam             |Alignments after QC and filtering. This is the set|
|                                       |used for subsequent analyses.                     |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.genome.bam               |Alignments of reads mapped on the chromosome      |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.trans.bam                |Alignments of reads reads mapped to transcripts   |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.mismapped.bam            |Alignments flagged as :term:`mismapped`           |
+---------------------------------------+--------------------------------------------------+

.. note::

   Please be aware that there is a difference between :term:`read` and :term:`alignment`
   counts. A single :term:`read` might map to several genomic positions and give rise
   to several alignments. Only if the pipeline employs unique-ness filtering will 
   :term:`read` and :term:`alignment` counts be the same.







