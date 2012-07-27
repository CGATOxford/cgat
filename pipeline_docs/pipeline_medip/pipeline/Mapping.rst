.. _Mapping:

=======
Mapping
=======

This section summarizes results of the read mapping process.

Overview
========

The following table summarizes the number of reads input, mapped and filtered.

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
   :slices: reads_total,reads_mapped,reverse,duplicate
   :tracks: r(.genome)

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reverse,duplicate
   :tracks: r(.genome)

   Mapping Results

Detailed results
================

Below are detailed results from the mapping process. It reports summary statistics from various
:term:`bam` formatted files that the pipeline creates.

.. toctree::
   :maxdepth: 2

   MappingAlignmentStatistics.rst
   CpGCoverage.rst

There are various alignment files created in the mapping process:

+---------------------------------------+--------------------------------------------------+
|*Filename*                             |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.genome.bam               |Alignments of reads mapped on the chromosome      |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.prep.bam                 |Alignments of reads after filtering               |
+---------------------------------------+--------------------------------------------------+

.. note::

   Please be aware that there is a difference between :term:`read` and :term:`alignment`
   counts. A single :term:`read` might map to several genomic positions and give rise
   to several alignments. Only if the pipeline employs unique-ness filtering will 
   :term:`read` and :term:`alignment` counts be the same.

