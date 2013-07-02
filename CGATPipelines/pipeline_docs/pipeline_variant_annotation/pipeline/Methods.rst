=======
Methods
=======

Gene sets
=========

ENESMBL
-------

The ensembl gene set is imported via the database dumps in :file:`.gtf` and :file:`.fa` format
taken from the *ftp* server.

Only protein-coding genes are considered. Note that there are no gene models on unplaced
contigs in these database dumps.

Refseq
------


Calling consequences
====================

Special cases
-------------

There are plenty of specials cases that the software needs to deal with.

1. ENSMUSG00000004099 (dnmt1): a stop codon appears in a transcript that is
   wrongly predicted. One transcript is longer than the others is likely to be
   over-predicted. Note that this stop codon is not part of a refseq transcript.

2. ENSMUSG00000022070: Indel at exon (99459458, 99459524) straddling intron - 
   recovers the splice site, but changes frame. Is it a true deletion?

3. The gene muc2 contains an X, becaus it starts with the sequence::

   >chr7:+:148930514:148930523                                                                                                                                                                                                                  
   NNAAGACCA


MGI
===

There are three options for importing MGI information:

1. Biomart 
    Downloading via biomart has the advantage, that also
    alleles without phenotypes can be obtained. However, there
    is a manual step involved. Includes alleles without phenotypes.
    The table is not normalized.

2. Database dumps
    Can be done automatically. However, the tables dumped are
    high-level views and not normalized. Also, not all alleles 
    are present, for example, ``MGI:4317524``. This information is
    in MRK_GeneTrap.rpt, but crucially, the MGI id missing.
    While the biomart table contains 282,383 alleles, I find
    only 38,411 in the database dumps. However, of the 282,383
    alleles, only 23,551 are alleles with phenotype.
    
    Note that some tables only contain only the high-level 
    phenotype.

3. SQL server
    Needs an account and requires specialized connectors to be 
    installed as their database is sybase. Also, to get meaningful
    reports the database schema has to be learned.

Prediction of deleterious SNPs
==============================

Polyphen
--------

We are using ``polyphen-2.0.18``.

+------------------------------+------------------------------+---------------------------------------------+
|parameter                     |value                         |comment                                      |
+------------------------------+------------------------------+---------------------------------------------+
|database                      |uniref100                     |                                             |
+------------------------------+------------------------------+---------------------------------------------+

The Polyphen pipeline was modified in the following way to use 
mouse sequences:

   * Instead of human annotations, mouse annotation was derived from swissprot/trembl entries
   * Instead of human annotations, mouse annotation was derived from PFAM
   * There is no dbSNP annotation.

The Polyphen classifier still uses the human training data. The developers suspect that this should
not make much of a difference.

Peptide sequences are first matched against the latest uniref100 release.
The matches are then scored.

In the `documentation <http://coot.embl.de/PolyPhen/>`_ it says::

   The user can choose between HumDiv- and HumVar-trained PolyPhen-2. Diagnostics of 
   Mendelian diseases requires distinguishing mutations with drastic effects from all 
   the remaining human variation, including abundant mildly deleterious alleles. Thus, 
   HumVar-trained PolyPhen-2 should be used for this task. In contrast, HumDiv-trained 
   PolyPhen-2 should be used for evaluating rare alleles at loci potentially involved in 
   complex phenotypes, dense mapping of regions identified by genome-wide association studies, 
   and analysis of natural selection from sequence data, where even mildly deleterious 
   alleles must be treated as damaging. 

We have used both ``HumVar`` and the ``HumDiv`` model.

In the analysis, SNPs are assigned a :term:`snp_id` and a :term:`locus_id`. 
The :term:`snp_id` refers to the SNP within a peptide sequence while the :term:`locus_id` refers
to the genomic location. 

If there are alternative transcripts overlapping a SNP, the same SNP will get two
:term:`snp_ids`, but the same :term:`locus_id`. As the peptide background might
be different for the same SNP depending on the transcript, its effect needs to be 
predicted twice. In either case, SNPs mapping to different bases will get separate
identifiers.

The distinction between a :term:`snp_id` and a :term:`locus_id` are important for
counting. For gene-based counts, use the :term:`locus_id`.

The document at `<http://sift.jcvi.org/www/annurev.genom.7.Ng_Henikoff.pdf>`_ contains
a review over various methods to predict deleterious SNPs. For Polyphen, it reports the 
false positive rate (FPR) to be 9% and the false negative rate (FNR) to be 31%. Note that 
the respective values for SIFT are FPR=20% and FNR=31%.

Panther
-------

We are using version ``csnpAnalysis1.01``.

The following parameters apply:

+------------------------------+------------------------------+---------------------------------------------+
|parameter                     |value                         |comment                                      |
+------------------------------+------------------------------+---------------------------------------------+
|database                      |PANTHER6.1                    |HMM database. Latest release as of 28.4.2010 |
+------------------------------+------------------------------+---------------------------------------------+
|substitution matrix           |blosum62                      |                                             |
+------------------------------+------------------------------+---------------------------------------------+
|dirichlet mixture             |9 component mixture           |                                             |
+------------------------------+------------------------------+---------------------------------------------+

Peptide sequences are aligned against the Panther HMM database using 
Panther default options (best hit ``-D b`` and use subfamilies ``-n``.

SNPs are then scored against the HMMs using the substitution matrix and
dirichlet prior listed in the table above.

Analysis of deleterious SNPs
----------------------------

Enrichment analysis
+++++++++++++++++++

The primary use case for variant prediction tools like polyphen_ and panther_
lies in the assessment of observed variants in one or few proteins. The question 
commonly asked is if a particular variant will affect or abrogate the protein's function.

In the genomic analysis of several strains, the question's focus has changed from the
variant to the gene set: Which genes are likely to have undergone functional change, 
given that we observe a collection of variants?

This question can be asked from two angles. The first angle investigates a single
strain with respect to the reference. The question is: which genes
are responsible for any of the phenotypic differences between strain and reference.
For this analysis, variants are collected and analysed per strain.

The second angle concerns the gene set itself. It asks, which genes are prone to
significant change. Here, variants are collected from all strains and analysed together.

In either case, the analysis needs to take into account that longer genes will contain
more SNPs simply due to random processes, and that there is a non-negligible false
positive rate for predicting deleterious variants. To avoid these effects, genes can
be flagged if they contain a larger proportion of deleterious SNPs than expected.
The expectation is given by a binomial distribution with a frequency corresponding
to the genomic SNP or deleterious SNP density overall.

The results of such an analysis need to be interpreted carefully. Genes might contain
an increased proportion of deleterious SNPs because of positive selection or relaxed
purifying selection. Similarly, the set might contain pseudogenes.

There are several enrichments that are relevant:
   * the proportion of SNPs in a gene compared to all codons in a gene (i.e., gene length)
   * the proportion of deleterious SNPs within SNPs
   * the proportion of deleterious SNPs in a gene compared to all codons in a gene:

Note that genes with only a single SNP which is deleterious are not likely to significant.
Similarly, short genes will not be called significant as there is not enough data.

Reference sequence
++++++++++++++++++

Is the fact that we are comparing against a reference sequence of importance?

Is there an ascertainment bias such that only variants to the reference sequence
are tested - but not the reference sequence positions? After all, they might be equally 
be considered a mutant if another strain was chosen as reference.

For example, I tested if there was an effect if I took a "probably damaging mutation"
from wildtype ``N`` to the mutant ``I`` and reverted the direction. The effects predicted were
identical after the feature collection:

+----------------+-------------------+------+---+---+---+---+------------------+--------------------+------+------+------+------+--------+--------+--------+--------+--------+                                                               
|snp_id          |     acc           |   pos|aa1|aa2|nt1|nt2|        prediction|            based_on|dScore|Score1|Score2|  Nobs| Nstruct|  IdPmax|  IdPSNP|  IdQmin|Comments|                                                               
+----------------+-------------------+------+---+---+---+---+------------------+--------------------+------+------+------+------+--------+--------+--------+--------+--------+                                                               
|snp0000213960   |ENSMUSP00000057998 |    90|  N|  I|  A|  T| probably damaging|           alignment|+2.775|+1.839|-0.936|   145|      20|   3.721|        |  77.419|        |                                                               
+----------------+-------------------+------+---+---+---+---+------------------+--------------------+------+------+------+------+--------+--------+--------+--------+--------+                                                               
|snp0000213960   |ENSMUSP00000057998r|    90|  I|  N|  T|  A| probably damaging|           alignment|-2.775|-0.936|+1.839|   145|      20|  31.143|  31.143|  89.677|        |                                                               
+----------------+-------------------+------+---+---+---+---+------------------+--------------------+------+------+------+------+--------+--------+--------+--------+--------+

Note that the PSIC scores are inverted, resulting in completely opposite delta score values. 

However, after running it through the classifier:

+----------------+-------------------+------+-----+-----+----------------+-------------------+------+---+---+---+---+------------------+--------------------+----------+------------------+----------+----------+----------+----------+--------+--------+--------+------+------+------+------+--------+------+------+-------+------+------+------+-------+------+------+------+------+------+--------+--------+--------+--------+--------+--------+--------+------+------+---+--------+------------+--------+--------+--------+--------+                                                                                                                                                                                               
|o_snp_id        |           o_acc   | o_pos|o_aa1|o_aa2|snp_id          |     acc           |   pos|aa1|aa2|nt1|nt2|        prediction|            based_on|    effect|        pph2_class| pph2_prob|  pph2_FPR|  pph2_TPR|  pph2_FDR|    site|  region|    PHAT|dScore|Score1|Score2|  Nobs| Nstruct| Nfilt|PDB_id|PDB_pos|PDB_ch| ident|length|NormAcc|SecStr|MapReg|  dVol| dProp|B-fact| H-bonds| AveNHet| MinDHet| AveNInt| MinDInt| AveNSit| MinDSit|Transv|CodPos|CpG| MinDJnc|     PfamHit|  IdPmax|  IdPSNP|  IdQmin|Comments|                                                                                                                                                                                               
+----------------+-------------------+------+-----+-----+----------------+-------------------+------+---+---+---+---+------------------+--------------------+----------+------------------+----------+----------+----------+----------+--------+--------+--------+------+------+------+------+--------+------+------+-------+------+------+------+-------+------+------+------+------+------+--------+--------+--------+--------+--------+--------+--------+------+------+---+--------+------------+--------+--------+--------+--------+                                                     
|snp0000213960   |ENSMUSP00000057998 |    90|    N|    I|snp0000213960   |ENSMUSP00000057998 |    90|  N|  I|  A|  T| probably damaging|           alignment|          |       deleterious|     0.973|    0.0613|     0.686|     0.151|        |        |        |+2.775|+1.839|-0.936|   145|      20|      |      |       |      |      |      |       |      |      |      |      |      |        |        |        |        |        |        |        |     1|     1|  0|        |            |   3.721|        |  77.419|        |                                                       
+----------------+-------------------+------+-----+-----+----------------+-------------------+------+---+---+---+---+------------------+--------------------+----------+------------------+----------+----------+----------+----------+--------+--------+--------+------+------+------+------+--------+------+------+-------+------+------+------+-------+------+------+------+------+------+--------+--------+--------+--------+--------+--------+--------+------+------+---+--------+------------+--------+--------+--------+--------+                                                                                                                                                                                               
|snp0000213960   |ENSMUSP00000057998r|    90|    I|    N|snp0000213960   |ENSMUSP00000057998r|    90|  I|  N|  T|  A|            benign|           alignment|          |           neutral|         0|     0.996|         1|     0.665|        |        |        |-2.775|-0.936|+1.839|   145|      20|      |      |       |      |      |      |       |      |      |      |      |      |        |        |        |        |        |        |        |     1|     1|  0|        |            |  31.143|  31.143|  89.677|        |                                                                                                                                                                                               
+----------------+-------------------+------+-----+-----+----------------+-------------------+------+---+---+---+---+------------------+--------------------+----------+------------------+----------+----------+----------+----------+--------+--------+--------+------+------+------+------+--------+------+------+-------+------+------+------+-------+------+------+------+------+------+--------+--------+--------+--------+--------+--------+--------+------+------+---+--------+------------+--------+--------+--------+--------+                                                                                                                                                                                               

The two substitutions are labeled appropriately. Thus, directionality in terms of the effect is built in.
However, the terms ``damaging`` and ``benign`` are relative to the choice of reference sequence. Again,
this makes sense in the original setting of Polyphen, where the sequence of a supposedly healthy individual
is contrasted with sequences of patients. It makes less sense if wild-type strains are compared.


Glossary
========

.. glossary::

   nmd
      `nonsense mediated decay <`http://en.wikipedia.org/wiki/Nonsense_mediated_decay>`_. Alleles
      are assumed to be affectd by nmd if they contain a stop-codon any exon but the last exon.

   variant
      A difference between the reference genome and the re-sequenced genome. A variant can
      be either a :term:`SNP` or an :term:`Indel`

   SNP
      A single nucleoutide polymorphism. Position at which the resequenced sequnce differs 
      from the reference sequence by the substitution of a single base.

   Indel
      Position at which the resequenced sequnce differs from the reference sequence by
      an insertion or deletion.

    dSNP
       A deleterious SNP. The deleteriousness is the result of a prediction.

   local effect
      A local effect of a :term:`variant` describes the consequence of this :term:`variant`
      within a transcript ignoring other variants in the same transcript. Local effects
      usually overestimate the effect of a variant. For example, a :term:`SNP` in position
      of 3 of a codon might create a stop codon, but a :term:`SNP` just two positions	
      upstream might correct the stop codon into a simple missense codon change.
