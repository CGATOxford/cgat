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
|reads_total                            |Number of reads submitted                         |
+---------------------------------------+--------------------------------------------------+
|reads_mapped                           |Number of reads mapped                            |
+---------------------------------------+--------------------------------------------------+
|mapped                                 |Number of alignments                              |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,mapped

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: matrix
   :slices: reads_total,reads_mapped,mapped
   :transform-matrix: normalized-row-max

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,mapped

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,mapped
   :transform-matrix: normalized-row-max

   Mapping Results

.. report:: Mapping.MappingFlagsHits
   :render: line-plot
   :as-lines:
   :layout: column-2
   :xrange: 0,10

   Number of reads per number of alignments (hits) per read.

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-2

   Number of mismatches per alignments.

FastQC
======

The following list presents links to the results of the :term:`fastqc` tool.
The fastqc tool is run over aligned reads only.

.. report:: Mapping.FastQCReport
   :render: user

   fastqc results







