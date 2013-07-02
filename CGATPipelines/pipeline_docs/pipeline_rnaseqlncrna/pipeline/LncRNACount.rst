========
LncRNA
========

Below are the descriptives for the lncRNA that have been predicted by the pipeline.


Number of LncRNA
==================

.. report:: LncRNACounts.Counts
   :render: table
   :groupby: all

   descriptives for predicted lncRNA and protein coding genes/transcripts


.. report:: LncRNACounts.Counts
   :render: interleaved-bar-plot
   :layout: column-2
   :groupby: all
   :slices: no_genes,no_transcripts,no_multi_exon_transcripts,no_single_exon_transcripts

   number of lncRNA and protein coding genes



Transcript and gene lengths
============================

.. report:: LncRNACounts.AveLength
   :render: table
   

   average length of lncRNA and protein coding genes


.. report:: LncRNACounts.AveLength
   :render: interleaved-bar-plot
   

   average length of lncRNA and protein coding transcripts/genes



Length distribution
====================

.. report:: LncRNACounts.LengthDistribution
   :render: line-plot
   :transform: histogram
   :as-lines:
   :layout: column-2


   length distribution for lncRNA and protein coding transcripts/genes
   

.. report:: LncRNACounts.LengthCumFreq
   :render: line-plot
   :as-lines:
   :layout: column-2
   
   length cumulative frequency for lncRNA and protein coding transcripts


Number of transcripts
======================

.. report:: LncRNACounts.TranscriptNumbers
   :render: table
   

   average number of transcripts at lncRNA and protein coding loci


.. report:: LncRNACounts.TranscriptNumbers
   :render: interleaved-bar-plot
   

   average number of transcripts at lncRNA and protein coding loci


.. report:: LncRNACounts.TranscriptNumberDistribution
   :render: histogram-plot
   :transform: histogram
   :groupby: track
   :as-lines:
   :layout: column-3

   transcript number distributions at lncRNA and protein coding loci


.. report:: LncRNACounts.TranscriptNumberCumFreq
   :render: line-plot
   :as-lines:

   transcript number cumulative frequencies for lncRNA and protein coding loci









