==========================
Number of allelic variants
==========================

This section shows the number of alleles per transcript.
There can be 0, 1 or two alleles depending wether both, one 
or none of the transcript alleles are identical to the reference
sequence.

.. report:: Effects.TranscriptsNumAlleles
   :render: stacked-bar-plot
   :slices: separate                                                                                                                                                                                                                         
   :layout: column-3
                                                                                                                                                                                                                                             
   Number of alleles per transcript.

.. report:: Effects.TranscriptsNumAlleles
   :render: table                                                                                                                                                                                                                            
   :slices: separate                                                                                                                                                                                                                         

   Number of alleles per transcript.

===================
Number of genotypes
===================

This section lists for each transcript whether the allelic variants
are resolvable. Resolvable variants are for example those in which all variants
are homozygous or all are heterozygous with one wild type allele. In the
latter case, for example, it is assumed that there is one wild-type and
one variant allele.

+------+-------------------+----------------------------------+----------+
| Code | Short explanation | Long explanation                 | Variants |
+------+-------------------+----------------------------------+----------+
|O     |Homozygous         |All variants are homozygous       |1         |
+------+-------------------+----------------------------------+----------+
|E     |Heterozygous       |Resolvable heterozygous, for      |2         |
|      |                   |example only one heterozygous     |          |
|      |                   |variant.                          |          |
+------+-------------------+----------------------------------+----------+
|W     |Heterozygous       |Resolvable heterozygous, with one |1         |
|      |                   |version being the wildtype.       |          |
+------+-------------------+----------------------------------+----------+
|V     |Variant            |Unresolvable heterozgygous, but   |1         |
|      |                   |wildtype allelle always present.  |          |
|      |                   |Resolved towards a single variant.|          |
+------+-------------------+----------------------------------+----------+
|v     |Variant            |Unresolvable heterozygous. Only   |2         |
|      |                   |one variant at each position, but |          |
|      |                   |a mixture of homosygous and       |          |
|      |                   |heterozygous positions.           |          |
+------+-------------------+----------------------------------+----------+
|M     |Mixture            |Mixture of heterozygous variants. |2         |
|      |                   |Can not be resolved. Variants are |          |
|      |                   |assigned randomly to alleles.     |          |
+------+-------------------+----------------------------------+----------+
|``-`` |Wildtype           |Transcript has no variants        |0         |
+------+-------------------+----------------------------------+----------+

.. report:: Effects.TranscriptsGenotypeResolvable
   :render: stacked-bar-plot
   :slices: separate                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Genotype of transcripts

.. report:: Effects.TranscriptsGenotypeResolvable                                                                                                                                                                                            
   :render: table                                                                                                                                                                                                                            
   :slices: separate                                                                                                                                                                                                                         
                          
   Genotype of transcripts                                                                                                                                                                                                                   

==============================
Number of transcripts per gene
==============================

The following plot simply presents the distribution of the number
of transcripts per gene.

.. report:: Effects.TranscriptsPerGene                                                                                                                                                                                                       
   :render: line-plot                                                                                                                                                                                                                          
   :transform: histogram                                                                                                                                                                                                                     
   :tf-range: 1,,1  
   :as-lines:

   Number of transcripts per gene

