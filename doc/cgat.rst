.. _cgat:

============================================
CGAT - Computational Genomics Analysis Tools
============================================

CGAT is a collection of tools for the computational genomicist written
in the python language. The tools have been developed and accumulated in various
genome projects (`Heger & Ponting, 2007`_, `Warren et al., 2008`_) and NGS pro-jects
(`Ramagopalan et al., 2010`_). The tools are continuously being developed
as part of the `CGAT Training programme`_.

The tools work from the command line, but can readily be installed
within frameworks such as `galaxy`_.

Please note that the tools are part of a larger code base also
including genomics and NGS pipelines. More information about these
is here (:ref:`here <contents>`).

Detailed instructions on installation, usage and tools are below,
followed by a quickstart guide.

.. toctree::
   :maxdepth: 2

   CGATInstallation
   CGATUsage
   CGATRecipes
   CGATReference

.. _quickstart:
   
Quickstart
==========

To install the CGAT tools, type::

   pip install CGAT

This will install the CGAT scripts and libraries together with the
required dependencies. See :ref:`Installation Instructions` for
dependencies and troubleshooting.

CGAT tools are run from the unix command line. Lets assume we have
the results of the binding locations of a ChIP-Seq experiment
(chipseq.hg19.bed) in bed format and we want to know, how many
binding locations are intronic, intergenic and within exons.

Thus, we need to create a set of genomic annotations denoting
intronic, intergenic regions, etc. with respect to a reference gene set.
Here, we download the GENCODE geneset (Harrow et al., 2012) in GTF
format from ENSEMBL (Flicek et al., 2013). 

The following unix statement downloads the ENSEMBL gene set containing
over-lapping transcripts, and outputs a set of non-overlapping genomic
annotations in gff format (:file:`annotations.gff`) by piping the data
through various GAT tools::
 
   wget -qO- ftp://ftp.ensembl.org/pub/release-72/gtf/homo_sapiens/Homo_sapiens.GRCh37.72.gtf.gz
   | gunzip
   | awk '$2 == "protein_coding"' 
   | gff2gff.py --genome-file=hg19 --sanitize=ucsc --skip-missing
   | gtf2gtf.py --sort=gene
   | gtf2gtf.py --merge-exons --with-utr
   | gtf2gtf.py --filter=longest-gene
   | gtf2gtf.py --sort=position
   | gtf2gff.py --genome-file=hg19 --flank=5000 --method=genome
   | gzip
   > annotations.gff.gz


.. note::
   The statements above need an indexed genome. To create such an
   indexed genome for hg19, type the following:
  
   wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
   | index_fasta.py hg19 - > hg19.log
   
CGAT tools can be chained into a single work flow using unix
pipes. The above sequence of commands in turn (1) reconciles UCSC and
ENSEMBL naming schemes for chromosome names, (2) merges all exons of
alternative transcripts per gene, (3) keeps the longest gene in case
of overlapping genes and (4) annotates exonic, intronic, intergenic
and flanking region (size=5kb) within and between genes.

Note that the creation of :file:`annotations.gff.gz` goes beyond
simple interval intersection, as gene structures have to be normalized
from multiple possible alternative transcripts to a single transcript
that is chosen by the user depending on what is most relevant for the
analysis.

Choosing different options can provide different sets of
answers. Instead of merging all exons per gene, the longest transcript
might be selected by replacing (2) with ``gtf2gtf.py
--filter=longest-transcript``. 
Or, instead of genomic annotations, regulatory domains such as defined by GREAT might be obtained by
removing (3) and replacing (4) with ``gtf2gff.py --method=great-domains``.

The generated annotations in annotations.gff can then be used to count
the number of transcription factor binding sites using bed-tools or
other interval intersections. Here, we will use another CGAT tool,
``gtf2table.py``, to do the counting and classification::

   zcat /ifs/devel/gat/tutorial/data/srf.hg19.bed 
   | bed2gff.py --as-gtf 
   | gtf2table.py --counter=classifier-chipseq --filename-gff=annotations.gff.gz

The scripts follow a consistent naming scheme centered around common
genomic formats. Because of the common genomic formats, the tools can
be easily combined with other tools such as bedtools (Quinlan and
Hall, 2010) or UCSC tools (Kuhn et al. 2013).


.. _Heger & Ponting, 2007: 
.. _Warren et al., 2008:
.. _Ramagopalan et al., 2010:
 
