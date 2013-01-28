.. _Mapping:

================
Mapping Overview
================

This section summarizes results of the mapping process.

Reads
=====

The following table summarizes the number of reads input, mapped and
filtered. Note that this table only looks at read numbers. For paired
end data, the number of reads are reported irrespective of whether
these were from the first or second pair.

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|reads_total                            |Number of reads submitted                         |
+---------------------------------------+--------------------------------------------------+
|reads_mapped                           |Number of reads mapped (after filtering)          |
+---------------------------------------+--------------------------------------------------+
|reads_norna                            |Number of reads not mapping to repetitive RNA     |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.MappingSummary
   :render: table
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping results in terms of reads

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping results in terms of reads

Paired reads
============

The following table summarizes the number of read pairs input and
mapped. Note that the number of pairs mapped can be larger than the 
number of pairs submitted as a pair might map to multiple locations.

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|pairs_total                            |Number of pairs submitted                         |
+---------------------------------------+--------------------------------------------------+
|pairs_mapped                           |Number of pairs mapped (proper pairs)             |
+---------------------------------------+--------------------------------------------------+

.. report:: Mapping.MappingSummary
   :render: table
   :slices: pairs_total,pairs_mapped

   Mapping results in terms of pairs

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: pairs_total,pairs_mapped

   Mapping results in terms of pairs





