.. _Mapping:

========
Overview
========

This section summarizes results of the read mapping process.

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

   Mapping Results

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments

   Mapping Results

.. note::

   Please be aware that there is a difference between :term:`read` and :term:`alignment`
   counts. A single :term:`read` might map to several genomic positions and give rise
   to several alignments. Only if the mapper employs unique-ness filtering will 
   :term:`read` and :term:`alignment` counts be the same.


Paired end data
===============

The following table summarizes the number of read pairs input and
mapped.

+---------------------------------------+--------------------------------------------------+
|*Column*                               |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|pairs_total                            |Number of reads submitted                         |
+---------------------------------------+--------------------------------------------------+
|reads_mapped                           |Number of reads mapped (after filtering)          |
+---------------------------------------+--------------------------------------------------+
|reads_norna                            |Number of reads not mapping to repetitive RNA     |
+---------------------------------------+--------------------------------------------------+






