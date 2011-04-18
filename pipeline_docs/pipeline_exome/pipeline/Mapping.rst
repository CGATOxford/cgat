=======
Mapping
=======

Mapping results
===============

+---------------------------------------+--------------------------------------------------+
|*Filename*                             |*Contents*                                        |
+---------------------------------------+--------------------------------------------------+
|:term:`track`.bam                      |Alignments after QC and filtering. This is the set|
|                                       |used for subsequent analyses.                     |
+---------------------------------------+--------------------------------------------------+


Alignments
----------

The following table present an overview of the alignments in the 
BAM files for each :term:`track`.

.. report:: Mapping.MappingSummary
   :render: table
   :slices: total,mapped,reverse,duplicates

   Mapping Summary

.. report:: Mapping.MappingSummary
   :render: interleaved-bar-plot
   :slices: total,mapped,reverse,duplicates

   Mapping Summary



