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
   :split-at: 10
   :layout: column-3

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

.. report:: Mapping.PairedMappingSummary
   :render: table
   :slices: pairs_total,pairs_mapped

   Mapping results in terms of pairs

.. report:: Mapping.PairedMappingSummary
   :render: interleaved-bar-plot
   :slices: pairs_total,pairs_mapped
   :split-at: 10
   :layout: column-3

   Mapping results in terms of pairs

Detailed summary
-----------------

The following table breaks down the number of reads mapped into several
categories.

+-------------------------+----------------------------------------+
|*Column*                 |*Content*                               |
+-------------------------+----------------------------------------+
|pairs_unmapped           |pairs in which neither read maps        |
+-------------------------+----------------------------------------+
|pairs_proper_unique      |pairs which are proper and map uniquely.|
+-------------------------+----------------------------------------+
|pairs_incomplete         |pairs in which one of the reads maps    |
|                         |uniquel, but the other does not map.    |
+-------------------------+----------------------------------------+
|pairs_proper_duplicate   |pairs which are proper and unique, but  |
|                         |marked as duplicates.                   |
+-------------------------+----------------------------------------+
|pairs_proper_multimapping|pairs which are proper, but map to      |
|                         |multiple locations.                     |
+-------------------------+----------------------------------------+
|pairs_not_proper_unique  |pairs mapping uniquely, but not flagged |
|                         |as proper                               |
+-------------------------+----------------------------------------+
|pairs_other              |pairs not in any of the above categories|
+-------------------------+----------------------------------------+

.. report:: Mapping.PairedMappingSummary
   :render: table
   :slices: pairs_proper_unique_alignments,pairs_unmapped,pairs_incomplete,pairs_not_proper_unique_alignments,pairs_other,pairs_proper_duplicate,pairs_proper_multimapping

   Mapping results in terms of pairs

.. report:: Mapping.PairedMappingSummary
   :render: stacked-bar-plot
   :slices: pairs_proper_unique_alignments,pairs_unmapped,pairs_incomplete,pairs_not_proper_unique_alignments,pairs_other,pairs_proper_duplicate,pairs_proper_multimapping
   :split-at: 10
   :layout: column-3

   Mapping results in terms of pairs





