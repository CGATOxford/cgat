===========
Mappability
===========

This section investigates the :term:`mappability` of various gene models.

The :term:`mappability` is estimated using :term:`mappability` tracks 
from the UCSC genome browse or elsewhere. Such tracks contain real-valued
data between 0 and 1 for each genomic position (or at some meaningful interval).
A 0 denotes that a sequence segment of certain length at a particular genomic 
position can not be mapped uniquely (for some level uniqueness), while a 1 
denotes that this position is uniquely mappable. Intermediate values may reflect
some level of degeneracy.

For human data, we are usually using the `CRG GEM alignability tracks
<http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeMapability>`_.

Genes and transcripts
=====================

The following plot shows the mean alignability for each transcrip, i.e, the average
mappability across all genomic positions within a transcript.

.. report:: Genemodels.GenesetMappability
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Distribution of mean mappability for transcipts in each geneset.

.. report:: Genemodels.GenesetMappability
   :render: box-plot

   Distribution of mean mappability for transcipts in each geneset.
