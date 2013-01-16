*****
OPTIC
*****

Purpose
-------

The optic pipeline performs orthology assignment in 
a group of species.

Setting up
----------

Download CGAT code
++++++++++++++++++

CGAT code can be obtained by checking the latest version out from
mercurial:

hg clone http://www.cgat.org/hg/cgat 

In the following, the location of the checked out code will be
referred to as <src>.

Requirements
++++++++++++

Optic requires the following tools to be installed.

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|postgres_           |                   |database                                        |
+--------------------+-------------------+------------------------------------------------+
|muscle_             |>=3.81.3           |Multiple alignment                              |
+--------------------+-------------------+------------------------------------------------+
|phyop_              |                   |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|treebest_           |>=0.1              |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|paml_               |>=4.4c             |evolutionary rate estimation                    |
+--------------------+-------------------+------------------------------------------------+

Genomes
++++++++

The first step is to build the files with the genomic
sequences. Download fasta files with genomic secquences
and index them using the <src>/IndexedFasta.py tools.

Store the genomes in a separate directory (referred to
later as <genome>).


Genesets
++++++++

The second step is to upload the genesets into the database. The
pipeline expects protein coding transcripts from ENSEMBL and requires
two files that can be downloaded from the `ENSEMBL ftp site <http://www.ensembl.org/info/data/ftp/index.html>`_:

* File with peptide sequences, usually called ``species``.pep.all.fa.gz
  
* File with exon coordinates in gtf format, usually called
    ``species``.gtf.gz

The following example shows how to upload the human gene set. We are
going to be using ensembl version 62 on hg19::

    mkdir -p genesets/hs62

    cd genesets/hs62

    ln -s <genomes>hg19.fasta genome.fasta

    ln -s <genomes>hg19.idx genome.idx

    ln -s <mirror>Homo_sapiens.GRCh37.62.gtf.gz reference.gtf.gz

    ln -s <mirror>Homo_sapiens.GRCh37.62.pep.all.fa.gz reference.pep.fa.gz

Create the makefile::

    python <src>setup.py -m ensembl -p cgat_hs62

Create database tables::

    make prepare

Upload data::

   make all

To verify all was ok, look at the file :file:`predictions.check`. This file
compares the ENSEMBL supplied peptide sequences with those that have
been uploaded into the database.

The data is now stored in the database schema cgat_hs62.

You need to do this for all species that you want to run OPTIC on.

Make sure that the naming is consistent with the genomes. Thus,
``hg19`` should both refer to the gene set for human, but also
to the genomic sequence files (:file:`hg19.fasta`) in the 
<genomes> directory.

Running OPTIC
-------------

The optic pipeline works from several directories.

Create a working directory::

    mkdir optic

Create a makefile::

    python <src>setup.py -m optic -p optic
    cd optic

Enter the :file:`data` directory and create the following files:

Makefile.inc
	common makefile file options. For example::

		## Global configuration options for OPTIC
		PARAM_PROJECT_NAME=cgat_proj007

		DIR_TMP=/tmp/

		PARAM_DIR_DATA=<optic>data/

		PARAM_SRC_SCHEMAS=cgat_hs62 cgat_mm62 cgat_gg62 cgat_ac62 cgat_xt62
		cgat_dr62 cgat_oa65

		PARAM_SPECIES_TREE=((((((cgat_hs62,cgat_mm62),cgat_oa65),cgat_gg62),cgat_ac62),cgat_xt62),cgat_dr62);

		PARAM_ANALYSIS_DUPLICATIONS_OUTGROUPS=cgat_dr62

		## CGAT cluster params

		DIR_SCRIPTS=<src>/

		PARAM_QUEUE=all.q
		PARAM_QUEUE_LOCAL=all.q
		PARAM_QUEUE_SERVER=all.q

schema2sp
	tsv file mapping species in database schema to swissprot name (required
	for treebest)::

	    # map of species names to swiss prot taxonomic names
	    # used for njtree
	    schema  sp
	    cgat_hs62       HUMAN
	    cgat_mm62       MOUSE
	    cgat_xt62       XENTR
	    cgat_dr62       DANRE
	    cgat_ac62       ANOCA
	    cgat_gg62       CHICK
	    cgat_oa65       ORNAN
    
species_tree
	the phylogeny of the species in newick format, for example::

	((((((cgat_hs62,cgat_mm62),cgat_oa65),cgat_gg62),cgat_ac62),cgat_xt62),cgat_dr62);
   
species_tree_permutations
	the phylogeny of the species and permutations of it. Usually
	just a copy of species_tree
	
files needed for web server:

schema2colour
	map species name to colour::

	    cgat_hs62       255,204,0
	    cgat_dr62       255,102,204
	    cgat_xt62       204,204,0
	    cgat_mm62       204,102,255
	    cgat_ac62       102,255,255
	    cgat_gg62       125,125,255
	    cgat_oa65       255,255,102

schema2name
	map species name to real name::

	    schema  name
	    cgat_hs62       H. sapiens
	    cgat_mm62       M. musculus
	    cgat_ac62       A. carolinensis
	    cgat_xt62       X. tropicalis
	    cgat_dr62       D. rerio
	    cgat_gg62       G. gallus
	    cgat_oa65       O. anatinus


schema2url
	map species name to ENSEMBL URL::

	    cgat_hs62	    displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_mm62       displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_ac62       displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_xt62       displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_dr62       displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_gg62       displayGene?schema=%(species)s&gene_id=%(gene)s
	    cgat_oa65       displayGene?schema=%(species)s&gene_id=%(gene)s

Preparing data
+++++++++++++++

Create export data files. This will create fasta file and exon
boundary files for all species::

   make -C export export_clustering

Phyop
+++++

Setup pairwise phyop runs and run them::

   make -C orthology_pairwise prepare

   nice -19 nohup make -C orthology_pairwise all

Wait a while...

Clustering
++++++++++

The clustering step combines pairwise orthology assignments into
clusters of potential orthologs::

   make -C orthology_multiple prepare
   make -C orthology_mulitple all


Multiple alignment
++++++++++++++++++

Next, multiple alignments are built for each cluster::

   make -C malis prepare
   make -C malis all
   make -C malis summary.dir
   make -C malis summary

Orthologous groups
+++++++++++++++++++

Based on the multiple alignments, trees are built within each cluster
and the trees are split into orthologous groups::


   make -C paralogy_trees prepare
   make -C paralogy_trees all
   make -C paralogy_trees analysis
   make -C paralogy_trees summary

Configuration
-------------

Edit the :file:`Makefile` to configure the pipeline.

.. _muscle: http://www.drive5.com/muscle/
.. _paml: http://abacus.gene.ucl.ac.uk/software/paml.html
.. _treebest: http://treesoft.sourceforge.net/treebest.shtml
.. _postgres: http://www.postgresql.org/
.. _phyop: http://www.cgat.org/~andreas/leotools-0.1.x86_64.linux.tar.gz
