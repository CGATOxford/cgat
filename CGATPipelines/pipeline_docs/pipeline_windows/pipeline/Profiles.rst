========================
Meta gene profiles
========================

Meta gene profiles are read densities computed across 
genomic features such as transcripts, transcription
start sites, exons, ...


Transcripts
===========

.. report:: Tracker.TrackerDataframes
   :render: ggplot
   :glob: transcriptprofiles.dir/*.separateexonprofilewithintrons.matrix.tsv.gz
   :aes: x='bin',y='background',color='track'
   :geom: geom_point(size=1)
   :split-at: 5
   :split-always: input
   :regex: .*/(.*).transcriptprofile.* 

   Background normalized read densities across meta-gene
   models, with exons separated into first, middle and last
   exons.

TSS profiles
============

.. report:: Tracker.TrackerDataframes
   :render: ggplot
   :glob: transcriptprofiles.dir/*.tssprofile.matrix.tsv.gz
   :aes: x='bin',y='background',color='track'
   :geom: geom_point(size=1)
   :split-at: 5
   :split-always: input

   Background normalized read densities around
   transcription start and end sites.

