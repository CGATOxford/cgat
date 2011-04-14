=======
Alleles
=======

This section uses the allele data set to predict the effect
of variants on transcripts. Knockouts can have effects on three
levels: the allele, the transcript and the gene. 

Alleles
-------

The following table counts the number of alleles that have been
affected due to a premature stop, a splice site disruption or
have been knocked out by :term:`nmd`.

+------------------------------+------------------------------+
|*name*                        |*content*                     |
+------------------------------+------------------------------+
|total                         |total number of alleles       |
+------------------------------+------------------------------+
|wildtype                      |alleles that are exactly like |
|                              |the wildtype                  |
+------------------------------+------------------------------+
|stop                          |alleles that are truncated due|
|                              |to premature stop codons      |
+------------------------------+------------------------------+
|nmd                           |alleles that have been knocked|
|                              |out due to :term:`nmd`        |
+------------------------------+------------------------------+
|splice                        |alleles that are truncated due|
|                              |to splice site disruption     |
+------------------------------+------------------------------+
|variants                      |alleles that are variants of  |
|                              |the transcript, but not       |
|                              |disrupted                     |
+------------------------------+------------------------------+

.. report:: Alleles.AlleleCounts                                                                                                                                                                                                 
   :render: table                                                                                                                                                                                                                            
   
   Counts per allele

Transcripts
-----------

The following table counts the number of transcripts that have been
affected due to stops, :term:`nmd` and splice site disruption. Note
that a transcript is disrupted only, if both of its alleles are affected.

+------------------------------+------------------------------+
|*name*                        |*content*                     |
+------------------------------+------------------------------+
|total                         |total number of transcripts   |
+------------------------------+------------------------------+
|wildtype                      |transrcipts that are exactly  |
|                              |like the wildtype (both       |
|                              |alleles)                      |
+------------------------------+------------------------------+
|truncated                     |transcripts that are truncated|
|                              |in both alleles or in which   |
|                              |one allele is truncated and   |
|                              |the other is knocked out      |
+------------------------------+------------------------------+
|knockout                      |transcripts that have been    |
|                              |knocked out in both alleles   |
+------------------------------+------------------------------+
|partially disrupted           |transcripts that have one     |
|                              |disrupted allele.             |
+------------------------------+------------------------------+
|variants                      |other variants                |
|                              |                              |
|                              |                              |
+------------------------------+------------------------------+

.. report:: Alleles.TranscriptsCounts                                    
   :render: table

   Counts per transcript		    

Genes
-----

The following table counts the number of genes that have been
affected due to stops, :term:`nmd` and splice site disruption. Note
that a gene is disrupted only, if all of its transcripts (and alleles) 
are affected. 

+------------------------------+------------------------------+
|*name*                        |*content*                     |
+------------------------------+------------------------------+
|total                         |total number of genes         |
+------------------------------+------------------------------+
|wildtype                      |genes that are exactly like   |
|                              |the wildtype (both alleles in |
|                              |all transrcipts)              |
+------------------------------+------------------------------+
|truncated                     |genes in which all alleles in |
|                              |all transcripts are truncated,|
|                              |or all transcripts are either |
|                              |truncated or knocked out.     |
|                              |                              |
+------------------------------+------------------------------+
|knockout                      |genes in which all alleles and|
|                              |transrcripts have been knocked|
|                              |out.                          |
+------------------------------+------------------------------+
|partially disrupted           |transcripts that have one     |
|                              |disrupted allele.             |
+------------------------------+------------------------------+
|variants                      |other variants                |
|                              |                              |
|                              |                              |
+------------------------------+------------------------------+

.. report:: Alleles.GenesCounts 
   :render: table 

   Counts per gene	    

.. note::
   The numbers in the three tables above do not sum up, because a wildtype
   allele is an allele without any difference to the wildtype, and not just	
   disruptive mutations.

Wildtypes
=========

The following table lists the number of wildtype alleles that
have been found per transcript

.. report:: Alleles.WildtypeCountsPerTranscript
   :render: pie-plot
   :layout: column-3
   :width: 300
   
   Wildtype alleles per transcript.

Knockouts
=========

NMD knockouts
-------------

The following table lists all the genes that have been knocked out due
to :term:`nmd`.

+------------------------------+---------------------------------+
|*name*                        |*content*                        |
+------------------------------+---------------------------------+
|gene_id                       |ENSEMBL gene id                  |
+------------------------------+---------------------------------+
|gene_name                     |gene name                        |
+------------------------------+---------------------------------+
|ntranscripts                  |number of transcripts            |
+------------------------------+---------------------------------+
|contig                        |contig of gene                   |
+------------------------------+---------------------------------+
|strand                        |strand of gene                   |
+------------------------------+---------------------------------+
|stops-start                   |start coordinates of stop codons |
+------------------------------+---------------------------------+
|stops-end                     |end coordinates of stop codons   |
+------------------------------+---------------------------------+

.. report:: Alleles.GenesNMDKnockouts 
   :render: table 

   Genes that have been knocked out.

The following table gives an overview of all the genes that have been
knocked out:

.. report:: Alleles.GenesNMDKnockoutsSummary
   :render: table 

   Genes that have been knocked out in any of the strains.


The following table gives an overview of all the genes that have been
knocked out together with their OMIM annotations:

.. report:: Alleles.GenesNMDKnockoutsWithOMIM
   :render: table
   :slices: omim_phenotype
   :force:         
   
   Genes that have been knocked out in any of the strains.

Single exon truncations
-----------------------

The following table lists all the single exon genes that have been knocked out.
A single exon gene is said to be knocked out if it is seriously truncated,
here, if no more that 40 amino acids remain of the translated sequence.

This is a simple criterion and could be improved. For example in the case of GPCRs a
gene could be said to have been knocked out if less than 7 transmembrane
segments remain.

+------------------------------+---------------------------------+
|*name*                        |*content*                        |
+------------------------------+---------------------------------+
|gene_id                       |ENSEMBL gene id                  |
+------------------------------+---------------------------------+
|gene_name                     |gene name                        |
+------------------------------+---------------------------------+
|ntranscripts                  |number of transcripts            |
+------------------------------+---------------------------------+
|cds_len                       |length of cds                    |
+------------------------------+---------------------------------+
|cds_original_l0en             |length of cds in wildtype        |
+------------------------------+---------------------------------+
|contig                        |contig of gene                   |
+------------------------------+---------------------------------+
|strand                        |strand of gene                   |
+------------------------------+---------------------------------+
|stops-start                   |start coordinates of stop codons |
+------------------------------+---------------------------------+
|stops-end                     |end coordinates of stop codons   |
+------------------------------+---------------------------------+

.. report:: Alleles.GenesSingleExonKnockouts
   :render: table

   Single exon genes that have been knocked out.


Splice truncations
------------------

The following table lists genes which have been truncated 
due to abrogation of splice sites:

+------------------------------+---------------------------------+
|*name*                        |*content*                        |
+------------------------------+---------------------------------+
|gene_id                       |ENSEMBL gene id                  |
+------------------------------+---------------------------------+
|gene_name                     |gene name                        |
+------------------------------+---------------------------------+
|ntranscripts                  |number of transcripts            |
+------------------------------+---------------------------------+
|contig                        |contig of gene                   |
+------------------------------+---------------------------------+
|strand                        |strand of gene                   |
+------------------------------+---------------------------------+

.. report:: Alleles.GenesSpliceTruncated
   :render: table 

   Genes that have been knocked out.

Note that the criteria applied are very strict, only genes that
are truncated in all transcripts are listed here.
