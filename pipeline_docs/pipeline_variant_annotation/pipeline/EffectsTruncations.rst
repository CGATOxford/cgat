=====================
Truncated transcripts
=====================

Nonsynonymous substitutions resulting in premature stop codons within coding 
regions can lead to 'knockout' of transcripts from the transcript set. 
This is because transcripts containing premature stop codons transcripts are likely to be degraded by
`nonsense mediated decay <http://en.wikipedia.org/wiki/Nonsense_mediated_decay>`_.

The minimum truncation is the smallest reduction in ORF length in both alleles. 
If the wildtype allele is present, the minimum truncation will be zero. 
The maximum truncation is the largest reduction in ORF length between two alleles. Even if the maximum truncation
is large, the transcript could still be functional, if the minimum truncation is small.

Premature stop codons are only considered to be effectual if they are more than 5
codons from the wild type stop codon.

Genes encoding selenoproteins have been removed.


List of truncated transcripts
+++++++++++++++++++++++++++++

The following table lists transcripts that have part of their coding sequence
removed due to premature stop codons. Only transcripts that have homozygous 
premature stop codons are included in the list below. 

+--------------------+--------------------------------------------+
|Column              |                                            |
+--------------------+--------------------------------------------+
|gene_id             |ENSEMBL gene id                             |
+--------------------+--------------------------------------------+
|transcript_id       |ENSEMBL transcript id                       |
+--------------------+--------------------------------------------+
|gene_name           |gene name                                   |
+--------------------+--------------------------------------------+
|cds_len             |length (in bases) of transcript             |
+--------------------+--------------------------------------------+
|truncated           |number of truncated codons                  |
+--------------------+--------------------------------------------+

.. report:: Knockouts.TranscriptListTruncatedStopsMin
   :render: table 
   :force: 
                                                     
   Table with transcripts that are truncated due to premature stop codons.
   Only homozygous variants are counted.

The following table lists transcripts that have part of their coding sequence
removed due to premature stop codons. The table includes all transcripts that have 
at least one allele with a premature stop codon.

.. report:: Knockouts.TranscriptListTruncatedStopsMax
   :render: table 
   :force:
                                                     
   Table with transcripts that are truncated due to premature stop codons.
   All variants are counted.


List of truncated genes
+++++++++++++++++++++++

The following table lists genes that have part of their coding sequence
removed due to premature stop codons. Only genes that have homozygous 
premature stop codons that affect *all* transcripts of a gene are included 
in the list below. 

+--------------------+--------------------------------------------+
|Column              |                                            |
+--------------------+--------------------------------------------+
|gene_id             |ENSEMBL gene id                             |
+--------------------+--------------------------------------------+
|gene_name           |gene name                                   |
+--------------------+--------------------------------------------+
|min(cds_len)        |minimum length (in bases) of any transcript |
+--------------------+--------------------------------------------+
|ntranscripts        |number of transcripts                       |
+--------------------+--------------------------------------------+
|truncated           |number of truncated codons                  |
+--------------------+--------------------------------------------+

.. report:: Knockouts.GeneListTruncatedStopsMin
   :render: table 
   :force: 
                                                     
   Table with genes that have truncation due to stops


