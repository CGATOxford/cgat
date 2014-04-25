.. _methods:

=======
Methods
=======

.. automodule:: pipeline_rnaseqdiffexpression
   :members:
   :inherited-members:
   :show-inheritance:

.. toctree::

   python/Trackers.rst


Data
====

The pipeline will pick up files ``*_export.txt.gz`` in the working directory. These files
are assumed to be in the solexa genome analyzer version 1 export format.  These will be converted 
via the samtools script :file:`export2sam.pl` into SAM format and converted/indexed into BAM files.


Combining replicates
====================

For combining replicates of several runs, the intervals between the two runs will be
intersected into now intervals. Attributes PeakVal,AvgVal,nProbes of the new intervals 
are re-computed from the original reads. Thus these values will differ from the 
original intervals that are based on binninng.

Expression data
===============

Differentially expressed genes
-------------------------------

Differentially expressed genes are identified using Welch's T-test
and a FDR-rate of 0.1.

Associating probesets with intervals
------------------------------------

There is no a-priori association of intervals from ChIP-Seq with probesets from
the expression data. We have to start from the assumption that intervals can
be anywhere in the vicinity of a gene and might inded overlap with the upstream
or possibly even downstream gene of its target. The exception might be intervals
that overlap the :term:`TSS` of a gene. 

With this in mind, I will create a graph an "m:n" association between probesets
and intervals not overlapping a :term:`TSS`, but will create a 1:1 association 
between probesets and intervals that do overlap a :term:`TSS`.

I will use a radius of 5kb on either side of a gene to make the association.

Note that the array has a preference for transcripts from refseq ids, while
I have a preference for ENSEMBL ids. In the transcripts track - both are available.
That should not be a problem as I want to link intervals to probeset ids.

Note also that there might be several probesets associated to a single gene.

Interesting threads about SAM analysis and the number of permutations:

https://stat.ethz.ch/pipermail/bioconductor/2007-November/020111.html

For p-values down to 0.001 you need 1000 permutations, which requires 6 samples.
(924 permutations?).

See :pmid:`18793455` on reproducibility. These authors recommend sorting
by fold and P-Value. If the number of replicates is small, the number of
permutations is limited, which introduces variability in the rank order
of P-Values alone.

Building a high-confidence set
==============================

The high confidence set contains genes that are
   * differentially expressed
   * have differentially occupied binding sites


Motif analysis
==============

MAST curves
-----------

MAST curves check the specificity of ChiP-Seq results against a known
motif. MAST estimates the probability of one or more motifs occuring in
each interval and returns an E-Value for each. The E-value indicates
the number of occurences of a motif in a database of random sequences
where each motif would match with a score at least as high as the highest
scoring motif in the sequence in question.

In MAST curves, the horizontal axis represents the E-Value. The green
curve *with_motifs* indicates the number of intervals that contain
a motif with that particular E-Value or less. The blue curve *explained*
is the difference of the curve *with_motifs* and the E-value and indicates
the number of intervals that can be explained by these intervals at a given
E-Value. The maximum of the *explained* curve 

See `Valouev et al. (2008) <http://www.nature.com/nmeth/journal/v5/n9/full/nmeth.1246.html>`_. 
MAST curves are explained in the `supplement <http://www.nature.com/nmeth/journal/v5/n9/extref/nmeth.1246-S1.pdf>`_.

See also `Kharchenko et al. (2008) <http://www.nature.com/nbt/journal/v26/n12/full/nbt.1508.html>`_
and their `supplement <http://www.nature.com/nbt/journal/v26/n12/extref/nbt.1508-S1.pdf>`_.

Significance of overlap
=======================

Significance of overlap
-----------------------

Parameters
==========

The following table lists all the pipeline parameters:

.. report:: Tracker.Config
   :tracks: pipeline.ini
   :render: table 

   Table with pipeline parameters.
