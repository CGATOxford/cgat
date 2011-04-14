==========
SNP counts
==========

Variant types
=============

Variant types are annotated with a one-letter code:

+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
| *code*        | *description*                                                        |
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|O              | homozygous substitution                                              |                                                                                                                                                                   
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|E              | heterozygous substitution                                            |
+---------------+----------------------------------------------------------------------+
|I              | insertion with respect to reference genome                           |
+---------------+----------------------------------------------------------------------+
|D              | deletion with respect to reference genome                            |
+---------------+----------------------------------------------------------------------+
|W              | wild type (reference genome)                                         |
+---------------+----------------------------------------------------------------------+

.. report:: Trackers.VariantTypeCounts
   :render: table

   Variant types in each set.

.. report:: Trackers.VariantTypeCounts
   :render: matrix
   :transform-matrix: normalized-row-total

   Variant types in each set.

.. report:: Trackers.VariantTypeCounts
   :render: pie-plot
   :groupby: track
   :layout: column-3
   :width: 300

   Variant types in each set.

Base annotation
===============

Bases in the genome are annotated according to their position with
respect to a reference gene set. Each base is assigned a code:

The codes are:

+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
| *code*        | *description*                                                        |
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|a              | first codon position within a complete codon                         |                                                                                                                                                                   
+---------------+----------------------------------------------------------------------+                                                                                                                                                                   
|b              | second codon position within a complete codon                        |
+---------------+----------------------------------------------------------------------+
|c              | third codon position within a complete codon                         |
+---------------+----------------------------------------------------------------------+
|d              | coding base, but in multiple frames or strands                       |
+---------------+----------------------------------------------------------------------+
|e              | non-coding base in exon                                              |
+---------------+----------------------------------------------------------------------+
|f              | frame-shifted base                                                   |
+---------------+----------------------------------------------------------------------+
|g              | intergenic base                                                      |
+---------------+----------------------------------------------------------------------+
|i              | intronic base                                                        |
+---------------+----------------------------------------------------------------------+
|l              | base in other RNA                                                    |
+---------------+----------------------------------------------------------------------+
|m              | base in miRNA                                                        |
+---------------+----------------------------------------------------------------------+
|n              | base in snRNA                                                        |
+---------------+----------------------------------------------------------------------+
|o              | base in snoRNA                                                       |
+---------------+----------------------------------------------------------------------+
|r              | base in rRNA (both genomic and mitochondrial)                        |
+---------------+----------------------------------------------------------------------+
|p              | base in pseudogene (including transcribed, unprocessed and processed)|
+---------------+----------------------------------------------------------------------+
|q              | base in retrotransposon                                              |
+---------------+----------------------------------------------------------------------+
|s              | base within a splice signal (GT/AG)                                  |
+---------------+----------------------------------------------------------------------+
|t              | base in tRNA (both genomic and mitochondrial)                        |
+---------------+----------------------------------------------------------------------+
|u              | base in 5' UTR                                                       |
+---------------+----------------------------------------------------------------------+
|v              | base in 3' UTR                                                       |
+---------------+----------------------------------------------------------------------+
|x              | ambiguous base with multiple functions.                              |
+---------------+----------------------------------------------------------------------+
|y              | unknown base                                                         |
+---------------+----------------------------------------------------------------------+


.. report:: Trackers.SNPBaseAnnotation
   :render: table

   Counts of SNPs in various positions

.. report:: Trackers.SNPBaseAnnotation
   :render: table
   :transform-matrix: normalized-row-total

   Percentage of SNPs in various positions

.. report:: Trackers.SNPBaseAnnotation
   :render: pie-plot
   :groupby: track
   :layout: column-3
   :width: 300

   Percentage of SNPs in various positions

Coding SNPs (those not lying within intronic or intergenic regions).
These will still include SNPs within non-coding transcripts, as long as
these are annotated in the reference gene set.

.. report:: Trackers.SNPBaseAnnotationWithoutNonCoding
   :render: matrix
   :transform-matrix: normalized-row-total

   Percentage of coding SNPs in various positions

.. report:: Trackers.SNPBaseAnnotationWithoutNonCoding
   :render: pie-plot
   :groupby: track
   :layout: column-3
   :width: 300

   Percentage of coding SNPs in various positions

Proportion of synonynymous versus non-synonymous SNPs
-----------------------------------------------------

