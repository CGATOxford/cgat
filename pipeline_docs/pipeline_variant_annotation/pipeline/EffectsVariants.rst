========================
Overview of variants
========================

This section looks at individual variant sites and
how they map to transcripts.

Number of variants per transcript
=================================

The following table lists the number of :term:`variants` 
per transcript within different strains. Note that some 
variants might be counted twice if they affect both splicing and the
coding sequence.

.. report:: Effects.VariantsPerTranscipt                                                         
   :render: line-plot                                                                            
   :transform: histogram                                                                         
   :tf-aggregate: normalized-total,cumulative                                                    
   :as-lines:                                                                                    
                                                                                                 
   Number of variants per transcript.

.. report:: Effects.VariantsPerTranscipt                                                         
   :render: table                                                                                
   :transform: stats                                                                             
                                                                                                 
   Number of variants per transcript.												 


CDS variants
============

Effect codes
------------

Effect codes per transcript
+++++++++++++++++++++++++++

Effect codes describe the :term:`local effect` of a variant. Counts
here are per transcript and snp. If a snp occurs in several transcript,
it will be counted several times.

+--------+-----------------------------------------------------+
|code    |description                                          | 
+--------+-----------------------------------------------------+
|F       |frameshifting effect                                 | 
+--------+-----------------------------------------------------+
|S       |synonymous effect                                    | 
+--------+-----------------------------------------------------+
|N       |nonsynonymous effect                                 | 
+--------+-----------------------------------------------------+
|I       |synonymous insertion - the inserted codon leaves     |
|        |the original codon intact                            |
+--------+-----------------------------------------------------+
|D       |synonymous deletion - the deleted codon leaves the   |
|        |original codon intact or removes it entirely         |
+--------+-----------------------------------------------------+
|X       |deleterious effect (stop codon)                      |
+--------+-----------------------------------------------------+

.. report:: Effects.VariantsCDSEffectCodes
   :render: table

   Number of local effects of variants

.. report:: Effects.VariantsCDSEffectCodes
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total   

   Number of local effects of variants

Effect codes per strain
+++++++++++++++++++++++++

Effect codes describe the :term:`local effect` of a variant. Counts
here are per snp and effect. If a SNP at a certain genomic position
has the same effect in several transcripts, it is only counted once.

.. report:: Effects.VariantsCDSEffectCodesPerStrain
   :render: table

   Number of local effects of variants

.. report:: Effects.VariantsCDSEffectCodesPerStrain
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total   

   Number of local effects of variants

Effect codes per position
+++++++++++++++++++++++++

Effect codes describe the :term:`local effect` of a variant. Counts
here are per snp and effect and strain. If a SNP at a certain genomic position
has the same effect in several transcripts in several strains, it is only counted 
once. Note that different effects at the same position might still be double 
counted.

.. report:: Effects.VariantsCDSEffectCodesPerPosition
   :render: table

   Number of local effects of variants

.. report:: Effects.VariantsCDSEffectCodesPerPosition
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total   

   Number of local effects of variants

Variant codes
-------------

Variant codes distinguish between SNPs and Indel variants.

+--------+-----------------------------------------------------+
|code    |description                                          |
+--------+-----------------------------------------------------+
|=       |single base substitution (SNP)                       |
+--------+-----------------------------------------------------+
|\+      |insertion variant (one or more bases)                |
+--------+-----------------------------------------------------+
|\-      |deletion variant (one or more bases)                 |
+--------+-----------------------------------------------------+

.. report:: Effects.VariantsCDSVariantCodes
   :render: table

   Number of substitution, insertion and deletion variants.

.. report:: Effects.VariantsCDSVariantCodes
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total   

   Number of substitution, insertion and deletion variants.

Variant types
-------------

Variant types describe whether a variant is homozygous, heterozygous, etc.

+--------+-----------------------------------------------------+
|code    |description                                          |
+--------+-----------------------------------------------------+
|E       |heterozygous SNP                                     |
+--------+-----------------------------------------------------+
|O       |homozygous SNP                                       |
+--------+-----------------------------------------------------+
|D       |homozygous deletion                                  |
+--------+-----------------------------------------------------+
|I       |homozygous insertion                                 |
+--------+-----------------------------------------------------+
|D,D     |heterozygous deletion, two non-wildtype alleles      |
+--------+-----------------------------------------------------+
|I,I     |heterozygous insertion, two non-wildtype alleles     |
+--------+-----------------------------------------------------+
|D,W     |heterozygous deletion, one wildtype allele           |
+--------+-----------------------------------------------------+
|I,W     |heterozygous insertion, one wildtype allele          |
+--------+-----------------------------------------------------+

.. report:: Effects.VariantsCDSVariantTypes
   :render: table

   Number of substitution, insertion and deletion variants.

.. report:: Effects.VariantsCDSVariantTypes
   :render: pie-plot
   :layout: column-4
   :width: 200

   Number of substitution, insertion and deletion variants.



