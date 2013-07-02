===============
Input overview
===============

This section summarizes results of the bam files that have been
submitted to the peak calling pipeline.

Overview
========

The following table summarizes the number of reads input, mapped and filtered.

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|reads_total                            |Total number of reads                             |
+---------------------------------------+--------------------------------------------------+
|reads_mapped                           |Number of reads mapped                            |
+---------------------------------------+--------------------------------------------------+
|duplicates                             |Number of duplicate alignments                    |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.BamSummary
   :render: table
   :slices: reads_total,reads_mapped,duplicates

   Mapping Results

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,duplicates

   Mapping Results

