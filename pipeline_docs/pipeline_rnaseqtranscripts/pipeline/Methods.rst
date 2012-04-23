.. _methods:

=======
Methods
=======

Classification of transcripts
=============================

Ab-initio prediction of transcript models yields a large number of transcripts.
To make sense of this very heterogeneous collection, transcripts are compared 
to the full reference gene set. The classification works as follows (see :file:`gtf2table.py`).

``Predicted`` transcripts are classified by checking overlap with all ``known`` transcripts in the reference
gene set.

For overlap of ``predicted`` transcript and ``known`` transcript a class is assigned:

+--------------------+------------------------------------------------------------+
|*Labels*            |*Contents*                                                  |
+--------------------+------------------------------------------------------------+
|complete            |Intron structure match - all exons are present, though first|
|                    |and last exon might be different.                           |
+--------------------+------------------------------------------------------------+
|extended-fragment   |At least one intron boundary shared. ``predicted``          |
|                    |transcript extends ``known`` transcript, but at the same is |
|                    |incomplete.                                                 |
+--------------------+------------------------------------------------------------+
|extension           |At least one intron boundary shared. ``predicted``          |
|                    |transcript extends ``known`` transcript.                    |
|                    |                                                            |
+--------------------+------------------------------------------------------------+
|fragment            |At least one intron boundary shared. ``predicted``          |
|                    |transcript is shorter than ``known`` transcript.            |
|                    |                                                            |
+--------------------+------------------------------------------------------------+
|alternative         |At least one intron boundary shared. ``predicted``          |
|                    |transcript has additional/missing exons/introns.            |
|                    |                                                            |
+--------------------+------------------------------------------------------------+
|unknown             |``predicted`` and ``known`` transcript overlap, but do not  |
|                    |fall in any of the above categories.                        |
|                    |                                                            |
+--------------------+------------------------------------------------------------+

Furthermore, ``predicted`` transcripts that do not overlap the exons of a ``known`` 
transcript, but only introns, are classed as ``intronic``. 

Transcripts that overlap neither introns nor exons of a ``known`` transcript are labeled ``intergenic``
if they are more than *5kb* away from a terminal exon. If they are within *5kb*, they are labeled 
``flank5`` or ``flank3`` depending on the orientation of the closest known transcript. If no strand
is present in either the transcript model nor the closest :term:`referenc` transcript, the class is ``flank``.

If the ``known`` transcript is protein coding, the ``predicted`` transcript is further 
checked if it overlaps with the ``known`` UTR only. These transcripts are labelled ``utr5``
and ``utr3``.

Additionally, the strandedness of the overlap is recorded as well (senes and antisense).

To decide, which transcript is the closest match, the resultant class are sorted as above 
in a priority list, where sense orientation has higher priority than anti-sense orientation.
For example, a ``predicted`` transcript will be rather labeled as ``sense intronic`` rather
than ``antisense fragment``.

The classification goes across overlapping genes thus giving preference to ``sense`` matches
compared to ``antisense`` matches.

Building a lincRNA set
======================

The :term:`lincRNA` set is a union of known :term:`lincRNA`s in the :term:`reference gene set`
and intergenic transcripts from the :term:`abinitio gene set` where transcripts in the :term:
`reference gene set` are given preference.

   * All transcripts from the :term:`reference gene set` that are annotated as ``lincRNA``
   * Transcripts from the :term:`abinitio gene set` that do not overlap (introns or exons)
       of transcripts annotated as ``protein_coding`` or ``lincRNA``.

Transcript models that overlap repetetive sequence are removed.

:term:`lincRNA` genes are often expressed at low level and thus the resultant transcript models 
are fragmentory. To avoid some double counting in downstream analyses, transcripts overlapping 
on the same strand are merged.

Transcripts need to have a length of at least 200 bp.

The lincRNA set is a mixture of several transcripts:

   * miscellaneous RNA, for example RMRP (RNA component of mitochondrial RNA processing endoribonuclease)

.. _CuffCompare:

Cuffcompare codes
=================

:term:`cuffcompare` compares the predicted transcripts from one or more experiments amongst each other. 
Similar transcripts are aggregated as :term:`transfrags` and designated a shared identifier permitting
similar transcripts to be tracked across experiments. Furthermore, each transfrag is compared to
the reference gene set and classified into the following categories:

+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|Priority|Code|Description                                                                                                                                              |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|1       | =  | Complete match of intron chain                                                                                                                          |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|2       | c  | Contained                                                                                                                                               |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|3       | j  | Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript                                                |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|4       | e  | Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron indicating a possible pre-mRNA fragment.                    |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|5       | i  | A transfrag falling entirely within a reference intron                                                                                                  |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|6       | o  | Generic exonic overlap with a reference transcript                                                                                                      |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|7       | p  | Possible polymerase run-on fragment (within 2Kbases of a reference transcript)                                                                          |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|8       | r  | Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|9       | u  | Unknown intergenic transcript                                                                                                                           |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|10      | x  | Exonic overlap with reference on the opposite strand                                                                                                    |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|11      | s  | An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)                                       |                              
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+                              
|12      | .  | indicates multiple classifications                                                                                                                      |
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+
|        | *  | Any of the above - used to display summaries on all transfrags                                                                                          |
+--------+----+---------------------------------------------------------------------------------------------------------------------------------------------------------+


Glossary
========

.. glossary::

   reference gene set
      The reference gene set. Usually the reference gene set is obtained from ENSEMBL_.
      The reference gene set contains protein-coding transcripts but also other genomic
      annotations such as :term:`lincRNA`, pseudogenes, small RNA, etc.

   abinitio gene set
      The ab-initio gene set is built from experimental data alone.

   lincRNA
      Long intergenic non-coding RNA.






