========================
Meta gene profiles
========================

Meta gene profiles are read densities computed across 
genomic features such as transcripts, transcription
start sites, exons, ...

Genes
=====

The following plots show the distribution of reads across the exonic
regions of protein coding genes (in 1000 bins) and 1kb of the
upstream/downstream region.

.. report:: Tracker.TrackerDataframes
   :render: ggplot
   :glob: transcriptprofiles.dir/*.geneprofile.matrix.tsv.gz
   :aes: x='bin',y='background',color='track'
   :geom: geom_point(size=1)
   :split-at: 5
   :split-always: input
   :groupby: all
   :regex: .*/(.*).transcriptprofile.* 

   Background normalized read densities across meta-gene
   models, with exons separated into first, middle and last
   exons.

With exons
==========

.. report:: Tracker.TrackerDataframes
   :render: ggplot
   :glob: transcriptprofiles.dir/*.separateexonprofilewithintrons.matrix.tsv.gz
   :aes: x='bin',y='background',color='track'
   :geom: geom_point(size=1)
   :split-at: 5
   :split-always: input
   :groupby: all
   :regex: .*/(.*).transcriptprofile.* 

   Background normalized read densities across meta-gene
   models, with exons separated into first, middle and last
   exons.

TSS profiles
============

.. report:: Tracker.TrackerDataframes
   :render: ggplot
   :glob: transcriptprofiles.dir/*.tssprofile.matrix.tsv.gz
   :regex: .*/(.*).transcriptprofile.*
   :aes: x='bin',y='background',color='track'
   :geom: geom_point(size=1)
   :groupby: all
   :split-at: 5
   :split-always: input

   Background normalized read densities around
   transcription start and end sites.

