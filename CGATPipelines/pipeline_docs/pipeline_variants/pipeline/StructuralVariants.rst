===================
Structural Variants
===================

Overview of structural variants
===============================

+--------------------+-----------------------------------------------+
|*Name*              |*Content*                                      |
+--------------------+-----------------------------------------------+
|event               |structural event (relative to reference)       |
+--------------------+-----------------------------------------------+
|name                |type of inserted/deleted sequence              |
+--------------------+-----------------------------------------------+
|counts              |numb r of events                               |
+--------------------+-----------------------------------------------+
|nval                |number of events with length information       |
+--------------------+-----------------------------------------------+
|mean                |mean number of bases of SVs                    |
+--------------------+-----------------------------------------------+
|median              |median number of bases of SVs                  |
+--------------------+-----------------------------------------------+
|sum                 |total number of bases of SVs                   |
+--------------------+-----------------------------------------------+

The following table shows the number of structural variants.

.. report:: StructuralVariants.StructuralVariantsSummary
   :render: table     
   :force: 
   :groupby: track

   Table with the number of structural variants in each strain.

SVs and Transcripts
===================

+--------------------+--------------------------+
|*Name*              |*Content*                 |
+--------------------+--------------------------+
|all                 |all structural variants   |
+--------------------+--------------------------+
|ins                 |only insertions           |
+--------------------+--------------------------+
|del                 |only deletions            |
+--------------------+--------------------------+
|other               |not insertions/deletions  |
+--------------------+--------------------------+

.. report:: StructuralVariants.StructuralVariantsSummaryTranscripts
   :render: table

   Number of transcripts overlapping various structural variants.

.. report:: StructuralVariants.StructuralVariantsSummaryTranscripts
   :render: interleaved-bar-plot

   Number of transcripts overlapping various structural variants.

SVs and genes
=============

+--------------------+-------------------------------------------------------------------+
|*Name*              |*Content*                                                          |
+--------------------+-------------------------------------------------------------------+
|naffected_ins       |Number of genes affected by an insertion                           |
+--------------------+-------------------------------------------------------------------+
|naffected_del       |Number of genes affected by a deletion                             |
+--------------------+-------------------------------------------------------------------+
|naffected_other     |Number of genes affected by other SV                               |
+--------------------+-------------------------------------------------------------------+
|naffected_all       |Number of genes affected by a SV                                   |
+--------------------+-------------------------------------------------------------------+
|ndeleted            |Number of genes with fully deleted transcripts                     |
+--------------------+-------------------------------------------------------------------+
|is_ins              |Number of genes with all transcripts affected by an insertion      |
+--------------------+-------------------------------------------------------------------+
|is_del              |Number of genes with all transcripts affected by a deletion        |
+--------------------+-------------------------------------------------------------------+
|is_other            |Number of genes with all transcripts affected by other SV          |
+--------------------+-------------------------------------------------------------------+
|is_all              |Number of genes with all transcripts affected by a SV              |
+--------------------+-------------------------------------------------------------------+
|is_deleted          |Number of genes with all transcripts fully deleted                 |
+--------------------+-------------------------------------------------------------------+

.. report:: StructuralVariants.StructuralVariantsSummaryGenes
   :render: table

   Table with number of genes affected by structural variants

Number of exons per gene that are fully deleted. Genes with a single exon
are likely candidates for reverse transposition.

.. report:: StructuralVariants.StructuralVariantsDeletedGenesExons
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-range: ,,1
   :yrange: 0,1
   :as-lines:

   Number of exons per gene that are fully deleted

List of genes that have been fully deleted.

.. report:: StructuralVariants.StructuralVariantsDeletedGenes
   :render: table
   :transform: group
   :tf-fields: gene_id

   List of genes that have been fully deleted.
