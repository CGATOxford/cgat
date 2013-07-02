=========
Knockouts
=========

Overview
========

Knockouts are genes/transcripts that are likely to be absent
from the transcript set due to premature stop codons or 
aberrant splicing. In the first section (:ref:`Truncations`),
such truncated transcripts and genes are listed. The
second section (:ref:`nonsense`) lists transcripts that
are likely to be removed due to 
`nonsense mediated decay <`http://en.wikipedia.org/wiki/Nonsense_mediated_decay>`_.

Genes encoding selenoproteins have been removed.

.. Truncations:

Truncated transcripts/genes 
============================

Truncation happens if a transcript contains a premature stop codon either
due to a SNP or due to an abrogated splice signal. The minimum truncation
is the minimum truncation within both alleles. 

It the wildtype is present, the truncation will be 0, unless the wild type contains a premature
stop codon. The maximum truncation is the maximum between two alleles. Even if the maximum truncation
is large, the transcript could still be functional, if the minimum truncation is small.

Premature stop codons are only considered to be effectual, if they are more than 
:param:`Knockouts.TranscriptsTruncatedStopsMin`
codons from the wild type codon.

Transcripts
-----------

The following tables/plots count truncations per transcript.
Truncations less than :param:`Knockouts.TranscriptsTruncatedStopsMin`

.. report:: Knockouts.TranscriptsTruncatedStopsMin                                                                                                                                                                                                   
   :render: table                                                                                                                                                                                                                            
   :slices: separate                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Number of transcripts with a minimum truncation due to stops.

.. report:: Knockouts.TranscriptsTruncatedStopsMin                                                                                                                                                                                                   
   :render: bar-plot                                                                                                                                                                                                                           
   :slices: separate                                                                                                                                                                                                                         
   :layout: column-3
                         																										                                                                                                                                                                                                                         
   Number of transcripts with a minimum truncation due to stops.

.. report:: Knockouts.TranscriptsTruncatedStopsMax                                                                                                                                                                                                 
   :render: table                                                                                                                                                                                                                            
   :slices: separate                                                                                                                                                                                                                         
                                                                                                                                                                                                                                             
   Number of transcripts with a maximum truncation due to stops.
   This shows the minimum truncation due to both alleles.

.. report:: Knockouts.TranscriptsTruncatedStopsMax                                                                                                                                                                                                 
   :render: bar-plot                                                                                                                                                                                                                           
   :slices: separate                                                                                                                                                                                                                         
   :layout: column-3
                                                                                                                                                                                                                                             
   Number of transcripts with a maximum truncation due to stops.
   This shows the maximum truncation due to any one allele.

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
                                                     
   Table with transcripts that are truncated due to premature stop codons.
   All variants are counted.

Genes
-----

A gene is truncated, if *all* of its transcripts have a minimum
truncation of at least :param:`Knockouts.GenesListTruncatedStopsMin`.

.. report:: Knockouts.GenesTruncatedStops                                                                                                                                                                                                      
   :render: table                                                                                                                                                                                                                            
                                                                                                                                                                                                                                             
   Number of genes with a minimun/maximum truncation due to stops

.. report:: Knockouts.GenesTruncatedStops                                                                                                                                                                                                      
   :render: interleaved-bar-plot
                                                                                                                                                                                                                                             
   Number of genes with a minimum/maximum truncation due to stops

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
|min(cds_len)        |minimum length (in bases )of any transcript |
|                    |                                            |
+--------------------+--------------------------------------------+
|ntranscripts        |number of transcripts                       |
+--------------------+--------------------------------------------+
|truncated           |number of truncated codons                  |
+--------------------+--------------------------------------------+

.. report:: Knockouts.GeneListTruncatedStopsMin
   :render: table 
   :force: 
                                                     
   Table with genes that have truncation due to stops

.. Nonsense:

Nonsense mediated decay
=======================

Transcipts
----------

The following tables/plots count the number of transcripts
that are likely to be affected by :term:`nmd`.

.. report:: Knockouts.TranscriptsNMDMin                   
   :render: table
   :slices: cds

   Number of transcripts knocked out by NMD

.. report:: Knockouts.TranscriptsNMDMin
   :render: bar-plot
   :slices: cds

   Number of transcripts knocked out by NMD

.. report:: Knockouts.TranscriptsNMDMax
   :render: table
   :slices: cds

   Number of transcripts affected by NMD

.. report:: Knockouts.TranscriptsNMDMax
   :render: bar-plot
   :slices: cds

   Number of transcripts affected by NMD

List of transcripts with NMD
+++++++++++++++++++++++++++++

The following table lists transcripts that are likely to be affected
by NMD.Only transcripts that have homozygous premature stop codons are 
included in the list below. Genes encoding selenoproteins have been removed 
as well.

+--------------------+--------------------------------------------+
|Column              |                                            |
+--------------------+--------------------------------------------+
|gene_id             |ENSEMBL gene id                             |
+--------------------+--------------------------------------------+
|transcript_id       |ENSEMBL transcript id                       |
+--------------------+--------------------------------------------+
|gene_name           |gene name                                   |
+--------------------+--------------------------------------------+
|min(cds_len)        |length (in bases) of transcript             |
+--------------------+--------------------------------------------+
|truncated           |number of truncated codons                  |
+--------------------+--------------------------------------------+
|last_exon_start     |start of last exon                          |
+--------------------+--------------------------------------------+

.. report:: Knockouts.TranscriptListNMDMin
   :render: table 
   :force: 
   :slices: cds
              
   Table with transcripts that are affected by NMD.

The following table lists transcripts that are likely to be affected
by NMD. This list includes all transcript, where at least one allele
is affected by :term:`nmd'. 

.. report:: Knockouts.TranscriptListNMDMax
   :render: table 
   :slices: cds
              
   Table with transcripts that have one allele that is affected by NMD.
