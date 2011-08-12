===============
Mapping results
===============

This section summarizes results of the read mapping process.

Overview
========

The following table summarizes the number of reads input, mapped and filtered.

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|mapped                                 |Number of reads mapped                            |
+---------------------------------------+--------------------------------------------------+
|duplicates                             |Number of duplicate alignments                    |
+---------------------------------------+--------------------------------------------------+
|reverse                                |Number of reads mapping to reverse strand         |
|                                       |                                                  |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.MappingSummary
   :render: table
   :slices: mapped,duplicates,reverse

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: mapped,duplicates,reverse

   Mapping Results

Detailed results
================

.. Below are detailed results from the mapping process. It reports summary statistics from various
.. :term:`bam` formatted files that the pipeline creates.

.. .. toctree::
..    :maxdepth: 2

..    MappingProcess.rst
..    MappingContext.rst
..    MappingAlignmentStatistics.rst

.. There are various alignment files created in the mapping process:

.. +---------------------------------------+--------------------------------------------------+
.. |*Filename*                             |*Contents*                                        |
.. +---------------------------------------+--------------------------------------------------+
.. |:term:`track`.bam                      |BAM files after mapping process.                  |
.. |                                       |                                                  |
.. +---------------------------------------+--------------------------------------------------+
.. |:term:`track`.prep.bam                 |BAM files after removing duplicate reads.         |
.. +---------------------------------------+--------------------------------------------------+
.. |:term:`track`.norm.bam                 |Normalized bam files                              |
.. +---------------------------------------+--------------------------------------------------+

.. .. note::

..    Please be aware that there is a difference between :term:`read` and :term:`alignment`
..    counts. A single :term:`read` might map to several genomic positions and give rise
..    to several alignments. Only if the pipeline employs unique-ness filtering will 
..    :term:`read` and :term:`alignment` counts be the same.

