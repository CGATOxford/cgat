===============
Mapping context
===============

This section summarizes the genomic context that reads have mapped to. This section
serves as a rough indicator to where reads have aligned too. The counting is naive:

* Counts are in terms of alignments. Thus a single read might contribute several counts
   to the same context or different contexts.

* Some genomic contexts can be overlapping, thus some alignments might be counted several
   times. 

* An alignment needs to map to a context over at least 50% of its bases. 
  Thus some alignments spanning several contexts might be dropped.

Context results
===============

The following table lists all genomic contexts that reads map to. 

.. report:: Mapping.MappingContext
   :render: table
   :force:

   Number of alignments that align in a certain genomic context

Ribosomal expression
====================

Ribosomal RNA is one of the most abundant transcripts in a cell and dominates RNASeq samples
until it is removed. The following plots and tables examine the number of alignments to
repetitive RNA. Repetetive RNA annotation is taken from the UCSC repeatmasker tracks.

.. report:: Mapping.MappingContext
   :tracks: r(.accepted)
   :render: table
   :slices: mapped,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA,ribosomal_coding

   Number of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track

.. report:: Mapping.MappingContext
   :tracks: r(.accepted)
   :render: pie-plot
   :pie-first-is-total: notRNA
   :groupby: track
   :slices: mapped,RNA,rRNA,scRNA,snRNA,srpRNA,tRNA,ribosomal_coding
   :layout: column-3
   :width: 200

   Proportion of alignments that align to repetitive RNA annotations from 
   the UCSC repeatmasker track

Protein coding expression
=========================

The following plots list the number of alignments to protein coding and (protein coding) 
pseudogene exons. The annotations are taken from the ENSEMBL gene set.

.. report:: Mapping.MappingContext
   :tracks: r(.accepted)
   :render: table
   :slices: mapped,protein_coding,pseudogene

   Number of alignments that align to protein coding genes or pseudo genes.

.. report:: Mapping.MappingContext
   :tracks: r(.accepted)
   :render: pie-plot
   :pie-first-is-total: genomic
   :groupby: track
   :slices: mapped,protein_coding,pseudogene
   :layout: column-3
   :width: 200

   Proportion of alignments that align to protein coding genes or pseudo genes.

