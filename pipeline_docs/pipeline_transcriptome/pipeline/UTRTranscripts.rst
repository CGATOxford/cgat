.. _UTRTranscripts:

======================
Non-coding transcripts
======================


UTR Transcripts
===============

UTR transcripts are exclusively located it the UTR of known
genes in the reference gene set. The following table shows
the number of transcripts within the UTR

.. report:: Trackers.KnownAnnotations
   :render: table
   :slices: known

   Annotations known transcripts.

.. report:: Trackers.UTRTranscripts
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportion of UTR transcripts in 5'/3' UTRs.

.. report:: Trackers.UTRTranscripts
   :render: matrix
   :transform-matrix: normalized-row-total

   Proportion of UTR transcripts in 5'/3' UTRs.

Direction of transcription
==========================

Checking the direction of transcription can give an indication
on the function of non-coding transcription.

Note that predicting the directionality of transcription might rely
on the presence of polyA tails or splice sites. These calls can be
erroneous. If neither polyA tails nor splice sites are available,
the direction might be determined by the directionality of reads,
which might bear no relationship to the direction of transcription.

Intronic transcripts
--------------------

.. report:: Trackers.TranscriptionDirectionIntronic
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportion of transcripts in introns that share or do
   not share the same direction of transcription compared to the 
   surrounding gene.

Intergenic transcripts
----------------------

.. report:: Trackers.TranscriptionDirectionIntergenic
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total

   Proportion of transcripts in intergenic sequences
   (within 10kb) that have the same or different direction 
   of transcription compared to the proximal gene.
