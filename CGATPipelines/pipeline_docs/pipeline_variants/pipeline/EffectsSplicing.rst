===================
Effects on splicing
===================

The following table contains information on transcripts that have variants.

+--------------------+-------------------------------------------------------------------------+
|Field               |Content                                                                  |
+--------------------+-------------------------------------------------------------------------+
|ncanonical          |number of canonical splice sites                                         |
+--------------------+-------------------------------------------------------------------------+
|nunchanged          |number of splice sites without variant                                   |
+--------------------+-------------------------------------------------------------------------+
|nnovel              |number of changes from a non-canonical splice site to a canonical one    |
+--------------------+-------------------------------------------------------------------------+
|nunknown            |number of changes between non-canonical splice sites                     |
+--------------------+-------------------------------------------------------------------------+
|nframeshifts        |number of frameshift introns                                             |
+--------------------+-------------------------------------------------------------------------+
|ncorrected_frames   |number of variants that correct frameshift introns                       |
+--------------------+-------------------------------------------------------------------------+
|nuncorrected_frames |number variants that do not correct frameshit introns                    |
+--------------------+-------------------------------------------------------------------------+
|nnonsynonymous      |number of variant changes between canonical splice sites                 |
+--------------------+-------------------------------------------------------------------------+
|synonymous          |number of variants that keep splice site intact                          |
+--------------------+-------------------------------------------------------------------------+
|ndisrupted          |number of variants that disrupt a splice site                            |
+--------------------+-------------------------------------------------------------------------+
|nintrons            |total number of splice sites                                             |
+--------------------+-------------------------------------------------------------------------+

.. report:: Effects.SplicingCounts                                                                                                                                                                                                           
   :render: table                                                                                                                                                                                                                            
                                                                                                                                                                                                                                             
   Number of splice sites and the variants affecting them.


Frame-shifting introns
======================

Frame-shifts are encoded as introns in the ENSEMBL gene set. The following
table list, how many have been affected and corrected by variants. 

Note that frame-shifting introns can be due both an in-frame stop-codon in the
reference sequence or a missing or inserted base or bases. In the first case,
the correction of a frame-shift might not be an indel, but a SNP changing the
codon away from the stop codon. 

.. report:: Effects.FrameShiftCorrection                                                                                                                                                                                                     
   :render: table                                                                                                                                                                                                                            
                                                                                                                                                                                                                                             
   Number of frameshifting introns and their variants. Counts
   are per variant. Thus a SNP appearing in multiple transcripts is
   counted multiple times.

The following plot shows the transcripts that have corrected frameshifts
in at least one strain.

.. report:: Effects.FrameShiftCorrectedTranscripts
   :render: matrix-plot
   :transform-matrix: transpose
   :palette: Blues

   Transcripts with corrected frameshifts.														    

