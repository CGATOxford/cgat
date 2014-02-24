===============
Mapping summary
===============

This section gives an overview over the mapping results.

Please be aware that there is a difference between :term:`read` and :term:`alignment`
counts. A single :term:`read` might map to several genomic positions and give rise
to several alignments. Only if the mapper employs unique-ness filtering will 
:term:`read` and :term:`alignment` counts be the same.

Alignment summary
=================

This section presents an overview of the :term:`alignments` in the 
BAM files for each :term:`track`. Note that a read can have multiple alignments.

.. report:: Mapping.BamSummary
   :render: table
   :slices: alignments_mapped,alignments_reverse,alignments_rna,alignments_duplicates

   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :slices: alignments_mapped,alignments_reverse,alignments_rna,alignments_duplicates
   :split-at: 10
   :layout: column-3
   
   Mapping summary

.. report:: Mapping.MappingFlagsMismatches
   :render: line-plot
   :as-lines:
   :layout: column-3
   :width: 200
   :split-at: 10

   Number of alignments per number of mismatches in alignment.

Reads summary
=============

This section presents an overview of the mapping results in terms 
of :term:`reads`.

.. report:: Mapping.BamSummary
   :render: table
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments
  
   Mapping summary

.. report:: Mapping.BamSummary
   :render: interleaved-bar-plot
   :slices: reads_total,reads_mapped,reads_norna,reads_norna_unique_alignments
   :split-at: 10
   :layout: column-3

   Mapping summary

.. report:: Mapping.MappingFlagsHits
   :render: line-plot
   :as-lines:
   :layout: column-3
   :width: 200
   :xrange: 0,10
   :split-at: 10

   Number of reads per number of alignments (hits) per read.

Fragment lengths
================

This section shows the fragment size distribution. 

.. report:: Mapping.PicardInsertSizeHistogram
   :render: line-plot
   :as-lines:
   :yrange: 0,
   :split-at: 10

   Histogram of fragment sizes

.. report:: Mapping.PicardInsertSizeMetrics
   :render: table
   :force:
   :split-at: 10

   Fragment size summary

Mapping qualities
=================

.. report:: Mapping.MappingQuality
   :render: line-plot
   :as-lines:
   :layout: column-3
   :width: 200
   :split-at: 10
   :logscale: y

   Distribution of mapping qualities
