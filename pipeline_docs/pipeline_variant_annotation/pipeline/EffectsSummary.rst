==========================================
Effect of Variants on protein-coding Genes 
==========================================

Effects per gene
++++++++++++++++

Effect codes describe the :term:`local effect` of a variant. Counts
here are per gene and effect. If a SNP at a certain genomic position
has the same effect in several transcripts, it is only counted once.

+--------+-----------------------------------------------------+
|Code    |Description                                          | 
+--------+-----------------------------------------------------+
|F       |Frameshift (indel)                                   | 
+--------+-----------------------------------------------------+
|S       |Synonymous substitution                              | 
+--------+-----------------------------------------------------+
|N       |Nonsynonymous substitution                           | 
+--------+-----------------------------------------------------+
|I       |synonymous Insertion - the inserted codon(s) leaves  |
|        |the original frame intact                            |
+--------+-----------------------------------------------------+
|D       |synonymous Deletion - the deleted codon(s) leaves    |
|        |the original frame intact                            |
+--------+-----------------------------------------------------+
|X       |Truncation (premature stop codon)                    |
+--------+-----------------------------------------------------+


.. report:: EffectsSummary.VariantsCDSEffectCodesPerStrain
   :render: table

   Number of local effects of variants

.. report:: EffectsSummary.VariantsCDSEffectCodesPerStrain
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Number of local effects of variants


Effects per transcript
++++++++++++++++++++++

Effect codes describe the :term:`local effect` of a variant. Counts
here are per transcript and snp. If a snp occurs in several transcript,
it will be counted several times.

+--------+-----------------------------------------------------+
|Code    |Description                                          | 
+--------+-----------------------------------------------------+
|F       |Frameshift (indel)                                   | 
+--------+-----------------------------------------------------+
|S       |Synonymous substitution                              | 
+--------+-----------------------------------------------------+
|N       |Nonsynonymous substitution                           | 
+--------+-----------------------------------------------------+
|I       |synonymous Insertion - the inserted codon leaves     |
|        |the original frame intact                            |
+--------+-----------------------------------------------------+
|D       |synonymous Deletion - the deleted codon leaves the   |
|        |original frame intact                                |
+--------+-----------------------------------------------------+
|X       |Truncation (premature stop codon)                    |
+--------+-----------------------------------------------------+

.. report:: EffectsSummary.VariantsCDSEffectCodes
   :render: table

   Number of local effects of variants

.. report:: EffectsSummary.VariantsCDSEffectCodes
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total   

   Number of local effects of variants




