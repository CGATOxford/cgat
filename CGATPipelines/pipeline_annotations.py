
"""===================
Annotation pipeline
===================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The annotation pipeline imports various third party annotations
or creates them for use in other pipelines.

The purpose of this pipeline is to automate and standardize the
way we retrieve and build genomic annotations but also to allow
sharing of annotations between projects and people. An important
part is the reconciliation of different data sources in terms
of chromosome names.

Common to all annotations in this pipeline is that they are genomic -
i.e. they are genomic intervals are relate to genomic intervals. Thus,
annotations are tied to a particular version of a genome. This pipeline
follows two principal releases: the UCSC_ genome assembly version and an
ENSEMBL_ geneset version.

The pipeline contains multiple sections that can be built on demand
or when relevant. Certain annotations (ENCODE, GWAS data) exist only
for specific species. The sections are:

assembly
   Genome assembly related information such as the location of
   gaps, chromosome lengths, etc.

ucsc
   Typical annotations downloaded from UCSC such as repeats.

ensembl
   The Ensembl gene set, reconciled with the assembly,
   and various subsets (coding genes, noncoding genes, ...).

geneset
   Annotations derived from the ENSEMBL gene set.

enrichment
   Annotations of genomic regions useful for enrichment
   analysis. These are derived from multiple input sources.

gwas
   GWAS data from the GWAS Catalog and DistlD

ontologies
   Ontology annotations (GO, KEGG) of genes.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The :file:`pipeline.ini` needs to be edited so that it points to the
appropriate locations of the auxiliary files. See especially:

1 section ``[ensembl]`` with the location of the ENSEMBL dump
    files (``filename_gtf``, filename_pep``, ``filename_cdna``)

2 section ``[general]`` with the location of the indexed genomic
    fasta files to use and the name of the genome (default=``hg19``),
    see :doc:`../modules/IndexedFasta`.

3 section ``[ucsc]`` with the name of the database to use (default=``hg19``).

Input
-----

This script requires no input within the :term:`working directory`, but
will look up some files in directories specified in the configuration
file :file:`pipeline.ini` and download annotations using mysql.

Running
-------

The pipeline can be run as any other CGAT pipeline, but as its purpose
is to provide a set of shared annotation between multiple projects
there is an etiquette to be followed:

Using the pipeline results
--------------------------

The annotations pipeline provides an interface for presenting its
results to other pipelines. The interface is defined in the file
:file:`pipeline.ini`. For example::

   [interface]
   # fasta file with cdna sequences
   cdna_fasta=ensembl.dir/cdna.fasta

The ini file of pipeline annotations can be loaded into the parameter
dictionary of your own pipeline::

    PARAMS.update(P.peekParameters(
         PARAMS["annotations_dir"],
         "pipeline_annotations.py",
         on_error_raise=__name__ == "__main__",
         prefix="annotations_"),
         update_interface=True)

Parameters from the annotation pipeline are now accessible via the
``annotations_`` prefix. As a result, the file
:file:`ensembl.dir/cdna.fasta` can be accessed as::

    PARAMS['annotations_cdna_fasta']

Extending the pipeline
-----------------------

Please feel free to add more annotations to the pipeline, but
considering its shared usage, please consult with others. In
particular, consider the following questions:

1. Is the annotation that I want to add genomic? For example,
   protein-protein interaction data should be organized separately.

2. Is the annotation of general interest? Do not add if an annotation
   is specific to a particular species or of very specialized
   interest. Note that there are some exceptions for annotations from
   certain species (human).

3. Is the annotation subjective? The pipeline consciously
   avoids providing annotations for regions such as promotors as their
   definition varies from person to person. Instead, the pipeline
   presents files with unambiguous coordinates such as transcription
   start sites. In the case of promotors, these could be derived from
   transcription start sites and the ``bedtools extend`` command.

4. What is the right format for the annotation? :term:`bed` formatted
   file are ideal for intervals with a single annotation. If multiple
   annotations are assigned with a feature, use :term:`gff`. For genes,
   use :term:`gtf`. Do not provide the same information with different
   formats - formats can be easily interconverted using CGAT tools.

Known problems
--------------

The pipeline takes its basic information about the genome and genes
from files downloaded from the genome browsers:

* UCSC: the genomic sequence in :term:`fasta` format.
* ENSEMBL: the gene set in :term:`GTF` format.

Additional data is downloaded from the genome browser databases either
via mysql or through biomart. It is thus important that the releases of
this additional data is consistent with the input files above.

.. note::

    The mechanism for getting biomart to download data for
    a particular ENSEMBL release involves changing the biomart server
    to an archive server.

Also, other data sources will have release cycles that are not tied
to a particular UCSC or ENSEMBL release. It is important to coordinate
and check when updating these other data sources.

Working with non-ENSEMBL species
--------------------------------

:doc:`pipeline_annotations` is very much wedded to annotations in ENSEMBL-
and UCSC_. Using a non-ENSEMBL species or non-UCSC species is possible by
building ENSEMBL- or UCSC-like input files. Even so, annotations that are
downloaded from the ENSEMBL or UCSC database will not be built. You will
thus need to ask if it is worth the effort.

As many other pipelines depend on the annotations in this pipeline it is
necessary to set up a :doc:`pipeline_annotations` stub. To do so, simply
build the config files by running::

   python <SRC>pipeline_annotations.py config

and create the files that are being used in the downstream pipeline
explicitely (for example, for protein coding genes)::

   mkdir ensembl.dir
   cp <MYDATADIR>/my_gtf_geneset.gtf.gz ensembl.dir/geneset_coding.gtf.gz

Roadmap
-------

There are many annotations that could possibly be brought into this pipeline:

* ENCODE data
     Can be used directly from a download directory?

* Genome segmentation based on ENCODE
     Definitions of enhancers, etc. Note that these will depend not on the
     genome, but on the cell type as well and thus might be project specific?

* Gene networks
     Functional assocation between genes. Outside of the
     scope of this pipeline?

* Mapability
     Mapability tracks are not available from all genomes. The pipeline
     could include runnig GEM on the assembly. For now it has been taken
     out as it is a rather long job.

Pipeline output
===============

The results of the computation are all stored in an sqlite relational
database file or as compressed files in genomic formats in the pipeline
directory. Output files are grouped by sections listed below.

The sections correspond to primary targets in the pipeline, i.e., to
build all annotations in the section ``assembly`` type::

   python <SRC>pipeline_annotations.py make assembly

Section: assembly
-----------------

Annotations derived from the genome assembly. Results are
in :file:`assembly.dir`.

contigs.tsv
   A :term:`tsv` formatted table with contig sizes

contigs.bed.gz
   bed file with contig sizes

contigs_ungapped.bed.gz
   :term:`bed` file with contigs excluding any gapped regions

gaps.bed.gz
   :term:`bed` file with gapped regions in contigs

genome.tsv.gz
   chromosome nucleotide composition and other stats

cpg.bed.gz
   filename with locations of CpG in bed format

gc_segmentation.bed.gz
   bed file with genome segmented into regions of similar G+C content
   using naive window based classification

gcprofile_bins.bed.gz
   bed file with genome segmented according to similar G+C content
   using the GCProfile method

Tables:
   genome
      Nucleotide composition of chromosomes

Section: ucsc
-------------

Various UCSC derived annotations. Results are in the
:file:`ucsc.dir`.

cpgislands.bed.gz
    :term:`bed` file with gapped regions in contigs

repeats.gff.gz
    :term:`gff` formatted file with structural/complex repeats

allrepeats.gff.gz
    :term:`gff` formatted file with all repeats including
    simple repeats

rna.gff.gz
    :term:`gff` formatted file with ribosomal rna annotations

repeats.gff.gz
    A :term:`gff` formatted file of repetitive sequences (obtained
    from UCSC repeatmasker tracks).

rna.gff.gz
    A :term:`gff` formatted file of repetitive RNA sequences in the genome
    (obtained from UCSC repeatmasker tracks).

mapability_xx.bed.gz
    Mapability files from UCSC CRG Alignability tracks. XX is the read
    length.

mapability_xx.bed.filtered.gz
    Similar to mapability_xx.bed.gz, but short regions of low mapability
    have been merged.

Tables
   repeat
       complex repeat locations
   repeat_counts
       Number of occurances for each repeat type.


Section: ensembl
----------------

Annotations within the ENSEMBL gene set after reconciliation
with the UCSC genome assembly. The results are in :file:`ensembl.dir`.
Annotations here are the original ENSEMBL annotations bar some
filtering.

geneset_all.gtf.gz
   The full gene set after reconciling with assembly. Chromosomes names are
   renamed to be consistent with the assembly and some chromosomes
   are optionally removed. This file is the starting point for
   all annotations derived from the ENSEMBL geneset.

geneset_cds.gtf.gz
   A :term:`gtf` formatted file with only the CDS parts of transcripts.
   This set will naturally include only coding transcripts. UTR regions
   have been removed.

geneset_exons.gtf.gz
   A :term:`gtf` formatted file with only the exon parts of transcripts.
   This set includes both coding and non-coding transcripts. Coding
   transcripts span both the UTR and the CDS.

geneset_coding_exons.gtf.gz
   :term:`gtf` file with exon parts of protein coding transcripts.
   All other features are removed. These are all features annotated
   as "protein_coding" in the ENSEMBL gtf file.

geneset_noncoding_exons.gtf.gz
   :term:`gtf` file with exon parts of non-coding transcripts
   all other features are removed. These are all transcripts not
   annotated as "protein_coding" in the ENSEMBL gtf file.

geneset_lincrna_exons.gtf.gz
   :term:`gtf` file with exon parts of lincRNA transcripts. These
   are transcripts annotated as "lincRNA" in the ENSEMBL gtf file.

geneset_flat.gtf.gz
   A :term:`gtf` formatted file of flattened gene
   models. All overlapping transcripts have been merged. This set
   includes both coding and non-coding transcripts.

peptides.fasta
   A :term:`fasta` formatted file of peptide sequences of coding
   transcripts.

cds.fasta
   A :term:`fasta` formatted file of coding sequence of coding transcripts.

cdna.fasta
   A :term:`fasta` formatted file of transcripts including both
   coding and non-coding parts.

Tables:
   transcript_info
       Information about transcripts (gene, biotype, status, ...)
       downloaded from biomart.

   transcript_synonyms
       Alternative names for transcripts

   gene_info
       Information about ENSEMBL genes in ENSEMBL gtf file.

   ensembl_to_entrez
       Table mapping ENSEMBL gene identifiers to
       ENTREZ identifiers

   cds_stats
       Table with nucleotide composition of each CDS in a transcript.

   gene_stats
        Table with nucleotide composition of each gene aggregated over
        all transcripts

   transcript_stats
        Table with nucleotide composition of each transcript.

   transcript_stats
        Table with amino acid composition of each protein product.

Section: geneset
----------------

Annotations derived from the ENSEMBL gene set. Annotations in
this section have been computed from the ENSEMBL gene set.
Results are in the directory :file:`geneset.dir`.

One group of :term:`bed` files outputs regions spanning whole
transcripts, genes, transcription start sites or transcription
termination sites for each of the gene sets build in
the ensembl section. These files are called
``<geneset>_<subset>_<region>.bed.gz``.

geneset
   coding, noncoding, lincrna

subset
   transcripts
      regions on a per-transcript level, multiple entries
      per gene, one for each transcript
   genes
      regions on a per-gene level, one entry for each gene

regions
   region
      the complete region spanning a transcript or gene
   tss
      the transcription start site a transcript. For genes,
      it is the most upstream TSS within a gene that is reported.
   tts
      the transcription termination site of a transcript. For genes,
      it is the most downstream TTS within a gene that is reported.
   tssregion
      the region spanning all TSS within a gene

The pipeline will also compute intergenic regions for each of these
datasets called ``<geneset>_intergenic.bed.gz``.

Note that the pipeline will not compute upstream or downstream flanks
or define promotor regions as the extend of these are usually project
specific. However, these files can be easily created using bed-tools
commands taking as input the files above.

Other files in this section are:

pseudogenes.gtf.gz
   A :term:`gtf` formatted file with pseudogenes. Pseudogenes are
   either taken from the ENSEMBL annotation or processed
   transcripts with similarity to protein coding sequence. As some
   protein coding genes contain processed transcripts without an
   ORF, Pseudogenes might overlap with protein coding transcripts
   This set is not guaranteed to be complete.

numts.gtf.gz
   set of potential numts. This set is not guaranteed to be complete.

Section: gwas
-------------

Data derived from GWAS databases. The files in this section represent
regions around SNPs that have been associated with certain traits or
diseases in GWAS experiments.

gwas_catalog.bed.gz
   :term:`bed` formatted file with intervals associated with various
   traits from the `gwas catalog`_. Regions are centered around the
   listed SNPs and extended by a certain amount.

gwas_distild.bed.gz
   :term:`bed` formatted file with LD blocks associated with various traits
   from the DistilD_ database


Note that the GWAS section is only available for human.

Section: ontologies
--------------------

Data in this section are ontology assignments for genes in the ENSEMBL
geneset.

go_ensembl.tsv.gz
   table with GO assignments for genes. GO assignments are downloaded
   from ENSEMBL.

goslim_ensembl.tsv.gz
   table with GOSlim assignments for genes

go_geneontology.tsv.gz
    table with terms from geneontology.org

go_geneontology_imputed.tsv.gz
    table with terms from geneontology.org, ancestral terms imputed.

kegg
    table with imported KEGG annnotations through biomart. Note
    that KEGG through this source might be out-of-date.

Tables:
   go_ensembl_assignments
      Table with GO assignments for each gene in ENSEMBL.

   goslim_ensembl_assignments
      Table with GOSlim assignments for each gene in ENSEMBL.

   kegg_assignments
      KEGG assignments

Section: enrichment
-------------------

This section contains useful files for genomic enrichment analysis
a la gat_. The annotations are derived from other annotations in
this pipeline. Output files are in the directory :file:`enrichment.dir`.

annotation_gff.gz
   A :term:`gff` formatted file annotating the genome with respect
   to the geneset.  Annotations are non-overlapping and are based
   only on protein coding transcripts.

genestructure.gff.gz
   A :term:`gff` file annotation genomic regions by gene structure

territories.gff.gz
   gff file with gene territories, .i.e. regions around protein
   coding genes.  Intergenic space between genes is split at the
   midpoint between two genes.

tssterritories.gff.gz
   gff file with tss territories

greatdomains.gff.gz
   gff file of regulator domains defined a la GREAT

genomic_context_bed=genomic_context.bed.gz
   bed-formatted file with genomic context

genomic_function_bed=genomic_function.bed.gz
   bed-formatted file with functional annotations

genomic_function_tsv=genomic_function.tsv.gz
   tsv-formatted file mapping terms to descriptions

Database design
---------------

Tables in the database usually represent genomic features such as
transcripts, genes or chromosomes. These are identified by the
following columns:

+--------------------+-----------------------------------------+
|*Column*            |*Content*                                |
+--------------------+-----------------------------------------+
|transcript_id       |ENSEMBL transcript identifier            |
+--------------------+-----------------------------------------+
|gene_id             |ENSEMBL gene id                          |
+--------------------+-----------------------------------------+
|contig              |Chromosome name                          |
+--------------------+-----------------------------------------+

For each :term:`bed`, :term:`gff` or :term:`gtf` file there is a
summary in the database called <file>_<format>_summary.  The summary
contains the number of intervals, nucleotides covered, etc. for that
particular file.

For :term:`gtf` files there is also a file with summary statistics
called <file>_gtf_stats.




Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz.

To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz
   tar -xvzf pipeline_annotations.tgz
   cd pipeline_annotations.dir
   python <srcdir>/pipeline_annotations.py make full

Code
====
"""
import sys
import shutil
import itertools
import csv
import re
import os
import collections
from ruffus import *
from bx.bbi.bigwig_file import BigWigFile
import sqlite3
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Database as Database
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineGO as PipelineGO
import CGATPipelines.PipelineBiomart as PipelineBiomart
import CGATPipelines.PipelineDatabase as PipelineDatabase
import CGATPipelines.PipelineUCSC as PipelineUCSC
import CGATPipelines.PipelineKEGG as PipelineKEGG
import CGAT.Intervals as Intervals


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


# Set parameter dictionary in auxilliary modules
PipelineGeneset.PARAMS = PARAMS
PipelineGO.PARAMS = PARAMS
PipelineDatabase.PARAMS = PARAMS
PipelineUCSC.PARAMS = PARAMS


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    return dbh


############################################################
# Assembly
@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_contigs'])
def buildContigSizes(infile, outfile):
    '''output contig sizes.
    '''
    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs = IOTools.openFile(outfile, "w")

    for contig, size in fasta.getContigSizes(with_synonyms=False).iteritems():
        outs.write("%s\t%i\n" % (contig, size))

    outs.close()


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_contigs_bed'])
def buildContigBed(infile, outfile):
    '''output bed file with contigs
    '''
    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs = IOTools.openFile(outfile, "w")

    for contig, size in fasta.getContigSizes(with_synonyms=False).iteritems():
        outs.write("%s\t%i\t%i\n" % (contig, 0, size))

    outs.close()


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       (PARAMS['interface_contigs_ungapped_bed'],
        PARAMS['interface_gaps_bed'],
        ))
def buildUngappedContigBed(infile, outfiles):
    '''output bed file with contigs.

    Genomic regions with gaps are excluded.
    '''
    prefix = P.snip(infile, ".fasta")
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs_nogap = IOTools.openFile(outfiles[0], "w")
    outs_gap = IOTools.openFile(outfiles[1], "w")
    min_gap_size = PARAMS["assembly_gaps_min_size"]

    for contig, size in fasta.getContigSizes(with_synonyms=False).iteritems():

        seq = fasta.getSequence(contig)

        def gapped_regions(seq):
            is_gap = seq[0] == "N"
            last = 0
            for x, c in enumerate(seq):
                if c == "N":
                    if not is_gap:
                        last = x
                        is_gap = True
                else:
                    if is_gap:
                        yield(last, x)
                        last = x
                        is_gap = False
            if is_gap:
                yield last, size

        last_end = 0
        for start, end in gapped_regions(seq):
            if end - start < min_gap_size:
                continue

            if last_end != 0:
                outs_nogap.write("%s\t%i\t%i\n" % (contig, last_end, start))
            outs_gap.write("%s\t%i\t%i\n" % (contig, start, end))
            last_end = end

        if last_end < size:
            outs_nogap.write("%s\t%i\t%i\n" % (contig, last_end, size))

    outs_nogap.close()
    outs_gap.close()


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS["interface_genome_tsv"])
def buildGenomeInformation(infile, outfile):
    '''compute genome composition information.'''

    statement = '''
    cat %(infile)s
    | python %(scriptsdir)s/fasta2table.py
        --section=length
        --section=cpg
    | gzip
    > %(outfile)s
    '''
    P.run()


@transform(buildGenomeInformation, suffix(".tsv.gz"), ".load")
def loadGenomeInformation(infile, outfile):
    '''load genome information.'''
    P.load(infile, outfile)

##################################################################
##################################################################
##################################################################
# build G+C segmentation
##################################################################


@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"]) + ".fasta",
       PARAMS["interface_gc_segmentation_bed"])
def buildGenomeGCSegmentation(infile, outfile):
    '''segment the genome into windows according to G+C content.'''

    statement = '''
    python %(scriptsdir)s/fasta2bed.py
        --method=fixed-width-windows-gc
        --window-size=%(segmentation_window_size)i
        --log=%(outfile)s.log
    < %(infile)s
    | python %(scriptsdir)s/bed2bed.py
        --method=bins
        --num-bins=%(segmentation_num_bins)s
        --binning-method=%(segmentation_method)s
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s'''

    P.run()


@follows(mkdir('assembly.dir'))
@files(os.path.join(PARAMS["genome_dir"],
                    PARAMS["genome"]) + ".fasta",
       "assembly.dir/gcprofile.bed.gz")
def runGenomeGCProfile(infile, outfile):
    '''segment the genome into windows according to G+C content.'''

    # on some cgat109 I got libstc++ error:
    # error while loading shared libraries: libstdc++.so.5
    # cannot open shared object file: No such file or directory
    to_cluster = False

    statement = '''
    cat %(infile)s
    | python %(scriptsdir)s/fasta2bed.py
        --verbose=2
        --method=GCProfile
        --gcprofile-min-length=%(segmentation_min_length)i
        --gcprofile-halting-parameter=%(segmentation_halting_parameter)i
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s
    '''
    P.run()


@merge(runGenomeGCProfile, PARAMS["interface_gc_profile_bed"])
def buildGenomeGCProfile(infile, outfile):
    '''aggregate windows with similar G+C content into bins.
    '''
    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/bed2bed.py
        --method=bins
        --num-bins=%(segmentation_num_bins)s
        --binning-method=%(segmentation_method)s
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s'''
    P.run()


@files(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"),
       PARAMS['interface_cpg_bed'])
def buildCpGBed(infile, outfile):
    '''build bed file with CpG locations.'''

    statement = '''
    python %(scriptsdir)s/fasta2bed.py
        --method=cpg
        --log=%(outfile)s.log
    < %(infile)s
    | bgzip
    > %(outfile)s
    '''

    P.run()

    statement = '''
    tabix -p bed %(outfile)s
    '''
    P.run()


# -----------------------------------------------------------------
# ENSEMBL gene set
@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"], PARAMS['interface_geneset_all_gtf'])
def buildGeneSet(infile, outfile):
    '''build a gene set - firstly, reconciles chromosome names by removing
       those that do not occur in the specified genome assembly;
       secondly, removes chromosome names specified in pipeline.ini
    '''

    statement = ['''zcat %(infile)s
    | grep 'transcript_id'
    | python %(scriptsdir)s/gff2gff.py
    --sanitize=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log ''']

    if PARAMS["ensembl_remove_contigs"]:
        # in quotation marks to avoid confusion with shell special
        # characters such as ( and |
        statement.append(
            ''' --remove-contigs="%(ensembl_remove_contigs)s" ''')

    statement.append(''' | gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()


@files(buildGeneSet, PARAMS['interface_geneset_flat_gtf'])
def buildFlatGeneSet(infile, outfile):
    '''build a flattened gene set.

    All transcripts in a gene are merged into a single transcript
    merging overlapping exons into a single exon.

    *infile* is an ENSEMBL gtf file.
    '''
    PipelineGeneset.buildFlatGeneSet(infile, outfile)


@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"], "ensembl.dir/gene_info.load")
def loadGeneInformation(infile, outfile):
    '''load the transcript set.'''
    PipelineGeneset.loadGeneInformation(infile, outfile)


@follows(mkdir('ensembl.dir'))
@files(buildFlatGeneSet, "ensembl.dir/gene_stats.load")
def loadGeneStats(infile, outfile):
    '''load the transcript set.'''
    PipelineGeneset.loadGeneStats(infile, outfile)


@files(buildGeneSet,
       PARAMS["interface_geneset_cds_gtf"])
def buildCDSTranscripts(infile, outfile):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS parts of exons are output - UTR's are removed.
    '''
    PipelineGeneset.buildCDS(infile, outfile)


@files(buildGeneSet,
       PARAMS["interface_geneset_exons_gtf"])
def buildExonTranscripts(infile, outfile):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildExons(infile, outfile)


@files(buildExonTranscripts, "ensembl.dir/transcript_stats.load")
def loadExonStats(infile, outfile):
    '''load the transcript set stats.'''
    PipelineGeneset.loadTranscriptStats(infile, outfile)


@files(buildGeneSet,
       PARAMS["interface_geneset_coding_exons_gtf"])
def buildCodingExonTranscripts(infile, outfile):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildCodingExons(infile, outfile)


@files(buildGeneSet,
       PARAMS["interface_geneset_noncoding_exons_gtf"])
def buildNonCodingExonTranscripts(infile, outfile):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildNonCodingExons(infile, outfile)


@files(buildGeneSet,
       PARAMS["interface_geneset_lincrna_exons_gtf"])
def buildLincRNAExonTranscripts(infile, outfile):
    '''build a collection of transcripts from the lincRNA
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildLincRNAExons(infile, outfile)


@transform((buildCDSTranscripts,
            buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           suffix(".gtf.gz"), "_gtf.load")
def loadTranscripts(infile, outfile):
    '''load the transcript set.'''
    PipelineGeneset.loadTranscripts(infile, outfile)


@files(buildCDSTranscripts, "ensembl.dir/cds_stats.load")
def loadCDSStats(infile, outfile):
    '''load the transcript set stats.'''

    PipelineGeneset.loadTranscriptStats(infile, outfile)


@follows(mkdir('ensembl.dir'))
@files(buildGeneSet, "ensembl.dir/transcript_info.load")
def downloadTranscriptInformation(infile, outfile):
    '''load transcript information.'''

    tablename = P.toTable(outfile)

    columns = {
        "ensembl_gene_id": "gene_id",
        "ensembl_transcript_id": "transcript_id",
        "ensembl_peptide_id": "protein_id",
        "gene_biotype": "gene_biotype",
        "transcript_biotype": "transcript_biotype",
        "source": "source",
        "status": "gene_status",
        "transcript_status": "transcript_status",
        "external_gene_id": "gene_name",
        "external_transcript_id": "transcript_name",
        "uniprot_sptrembl": "uniprot_id",
        "uniprot_genename": "uniprot_name",
    }

    # biomart db dmelanogaster_gene_ensemble doesn't have attribute
    # uniprot_genename
    if PARAMS["genome"].startswith("dm"):
        del columns["uniprot_genename"]
    # same fix for yeast
    if PARAMS["genome"].startswith("sac"):
        del columns["uniprot_genename"]

    data = PipelineBiomart.biomart_iterator(
        columns.keys(),
        biomart=PARAMS[
            "ensembl_biomart_mart"],
        dataset=PARAMS[
            "ensembl_biomart_dataset"],
        host=PARAMS["ensembl_biomart_host"])

    # The full list of genes from this table is too extensive. The
    # following are removed: 1. Some genes are present as LRGxxx
    # identifiers """LRG stands for Locus Reference Genomic. An LRG is
    # a fixed sequence, independent of the genome, specifically
    # created for the diagnostic community to record DNA sequence
    # variation on a fixed framework""" These are removed below:
    data = filter(lambda x: not x['ensembl_gene_id'].startswith("LRG"), data)

    # 2. Some genes are present on artificial chromosomes such as
    # ENSG00000265928 on HG271_PATCH.
    # To filter these out, the gene ids are cross-checked against those in
    # the ensembl gtf file.
    gene_ids = set()
    with IOTools.openFile(infile) as inf:
        for gtf in GTF.iterator(inf):
            gene_ids.add(gtf.gene_id)

    data = filter(lambda x: x['ensembl_gene_id'] in gene_ids, data)

    PipelineDatabase.importFromIterator(
        outfile, tablename,
        data,
        columns=columns,
        indices=("gene_id", "transcript_id", "protein_id",
                 "gene_name", "transcript_name", "uniprot_id")
    )

    # validate: 1:1 mapping between gene_ids and gene_names
    dbh = connect()

    data = Database.executewait(dbh, """
    SELECT gene_name, count(distinct gene_id) from %(tablename)s
    GROUP BY gene_name
    HAVING count(distinct gene_id) > 1""" % locals())

    l = data.fetchall()
    if len(l) > 0:
        E.warn("there are %i gene_names mapped to different gene_ids" % len(l))
    for gene_name, counts in l:
        E.info("ambiguous mapping: %s->%i" % (gene_name, counts))

    # adding final column back into transcript_info for dmelanogaster genomes
    if PARAMS["genome"].startswith("dm"):
        Database.executewait(
            dbh,
            '''ALTER TABLE Table1 ADD COLUMN uniprot_name NULL''')

    # adding final column back into transcript_info for scerevisiae genomes
    if PARAMS["genome"].startswith("sac"):
        Database.executewait(
            dbh,
            '''ALTER TABLE transcript_info ADD COLUMN uniprot_name NULL''')


@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"],
       "ensembl.dir/ensembl_to_entrez.load")
def downloadEntrezToEnsembl(infile, outfile):
    '''load table to convert from ENSEMBL gene ids to entrez gene ids'''

    tablename = P.toTable(outfile)

    columns = {
        "ensembl_gene_id": "gene_id",
        "entrezgene": "entrez_id"}

    data = PipelineBiomart.biomart_iterator(
        columns.keys(),
        biomart="ensembl",
        dataset=PARAMS["ensembl_biomart_dataset"])

    PipelineDatabase.importFromIterator(
        outfile,
        tablename,
        data,
        columns=columns,
        indices=("gene_id", "entrez_id"))


@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_gtf"],
       "ensembl.dir/transcript_synonyms.load")
def downloadTranscriptSynonyms(infile, outfile):
    '''load table with synonyms for transcript identifiers.'''

    tablename = P.toTable(outfile)

    columns = {
        "ensembl_transcript_id": "transcript_id",
        "external_transcript_id": "transcript_name",
        "refseq_mrna": "refseq_id",
    }

    data = PipelineBiomart.biomart_iterator(
        columns.keys(),
        biomart=PARAMS[
            "ensembl_biomart_mart"],
        dataset=PARAMS[
            "ensembl_biomart_dataset"],
        host=PARAMS["ensembl_biomart_host"])

    PipelineDatabase.importFromIterator(
        outfile,
        tablename,
        data,
        columns=columns,
        indices=(
            "transcript_id", "transcript_name", "refseq_id"))


@follows(mkdir('ensembl.dir'))
@files(((PARAMS["ensembl_filename_pep"],
         PARAMS["interface_peptides_fasta"]), ))
def buildPeptideFasta(infile, outfile):
    '''load ENSEMBL peptide file

    *infile* is an ENSEMBL .pep.all.fa.gz file.
    '''
    PipelineGeneset.buildPeptideFasta(infile, outfile)


@follows(mkdir('ensembl.dir'))
@files(((PARAMS["ensembl_filename_cdna"],
         PARAMS["interface_cdna_fasta"]), ))
def buildCDNAFasta(infile, outfile):
    '''load ENSEMBL peptide file

    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''
    PipelineGeneset.buildCDNAFasta(infile, outfile)


@follows(mkdir('ensembl.dir'))
@merge(buildCDSTranscripts, PARAMS["interface_cds_fasta"])
def buildCDSFasta(infile, outfile):
    '''build cds sequences from peptide and cds file.

    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''

    PipelineGeneset.buildCDSFasta(infile, outfile)


@follows(mkdir('ensembl.dir'))
@files(PARAMS["ensembl_filename_pep"],
       "ensembl.dir/protein_stats.load")
def loadProteinStats(infile, outfile):
    '''load the transcript set.'''

    PipelineGeneset.loadProteinStats(infile, outfile)


@merge((loadProteinStats, downloadTranscriptInformation),
       "ensembl.dir/seleno.list")
def buildSelenoList(infile, outfile):
    '''export a list of seleno cysteine transcripts.'''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''
    SELECT DISTINCT transcript_id
    FROM transcript_info as t,
         protein_stats as p
    WHERE p.protein_id = t.protein_id AND
         p.nU > 0
    '''
    outf = open(outfile, "w")
    outf.write("transcript_id\n")
    outf.write("\n".join(
        [x[0]
         for x in Database.executewait(dbh, statement)]) + "\n")
    outf.close()


# ---------------------------------------------------------------
# geneset derived annotations
@follows(mkdir('geneset.dir'))
@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_region.bed.gz')
def buildTranscriptRegions(infile, outfile):
    '''Bed file of regions spanning a whole transcript.'''
    statement = """
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_region.bed.gz')
def buildGeneRegions(infile, outfile):
    '''Bed file of regions spanning whole genes.'''
    statement = """
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_tss.bed.gz')
def buildTranscriptTSS(infile, outfile):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1. '''
    statement = """
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --join-exons
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py --method=promotors
    --promotor=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_transript_tts.bed.gz')
def buildTranscriptTTS(infile, outfile):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1. '''
    statement = """
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --join-exons
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py --method=tts
    --promotor=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tss.bed.gz')
def buildGeneTSS(infile, outfile):
    '''annotate transcription termination sites from reference gene set. '''
    statement = """gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s"""
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tts.bed.gz')
def buildGeneTTS(infile, outfile):
    '''annotate transcription termination sites from reference gene set. '''
    statement = """gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py --method=tts --promotor=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s"""
    P.run()


@transform((buildCodingExonTranscripts,
            buildNonCodingExonTranscripts,
            buildLincRNAExonTranscripts),
           regex('.*geneset_(.*)_exons.gtf.gz'),
           r'geneset.dir/\1_gene_tssinterval.bed.gz')
def buildGeneTSSInterval(infile, outfile):
    '''create a single interval that encompasses all annotated TSSs
    for a given gene
    '''
    statement = """
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1
    --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
    | sed s/transcript/exon/g
    | sed s/exon_id/transcript_id/g
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id
    --log=%(outfile)s.log
    | gzip
    > %(outfile)s """
    P.run()


@transform(buildGeneRegions,
           regex('(.*)_.*.bed.gz'),
           add_inputs(buildContigSizes),
           r'\1_intergenic.bed.gz')
def buildIntergenicRegions(infiles, outfile):
    ''' Genomic regions not associated with genes.
    '''

    infile, contigs = infiles

    statement = '''zcat %(infile)s
    | complementBed -i stdin -g %(contigs)s
    | gzip
    > %(outfile)s'''
    P.run()


# ---------------------------------------------------------------
# UCSC derived annotations
@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_rna_gff"]), ))
def importRNAAnnotationFromUCSC(infile, outfile):
    '''import repetetive RNA from a UCSC formatted file.

    The repetetive RNA are taken from the repeat-masker track.

    The results are stored as a :term:`gff` formatted file.
    '''

    repclasses = "','".join(PARAMS["ucsc_rnatypes"].split(","))
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC(dbhandle, repclasses, outfile)


@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_repeats_gff"]), ))
def importRepeatsFromUCSC(infile, outfile):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses = P.asList(PARAMS["ucsc_repeattypes"])
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC(dbhandle, repclasses, outfile)


@follows(mkdir('ucsc.dir'))
@files(((None, PARAMS["interface_cpgislands_bed"]), ))
def importCpGIslandsFromUCSC(infile, outfile):
    '''import cpg islands from UCSC

    The repeats are stored as a :term:`bed` formatted file.
    '''

    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getCpGIslandsFromUCSC(dbhandle, outfile)


@transform(importRepeatsFromUCSC, suffix(".gff.gz"), ".gff.gz.load")
def loadRepeats(infile, outfile):
    '''load total repeats length'''

    headers = "contig,start,stop,class"
    statement = """zcat %(infile)s
    | python %(scriptsdir)s/gff2bed.py --name=class | grep -v "#"
    | cut -f1,2,3,4
    | python %(scriptsdir)s/csv2db.py
    --retry
    --table=repeats
    --header=%(headers)s
    --index=class
    > %(outfile)s; """
    P.run()


@transform(loadRepeats, suffix(".gff.gz.load"), ".counts.load")
def countTotalRepeatLength(infile, outfile):
    ''' Count total repeat length and add to database '''
    dbhandle = sqlite3.connect(PARAMS["database"])
    cc = dbhandle.cursor()
    statement = """DROP TABLE IF EXISTS repeat_length"""
    Database.executewait(dbhandle, statement)

    statement = """create table repeat_length as
    SELECT sum(stop-start) as total_repeat_length from repeats"""
    Database.executewait(dbhandle, statement)

    P.touch(outfile)


@follows(mkdir('uscs.dir'))
@files(((None, PARAMS["interface_allrepeats_gff"]), ))
def importAllRepeatsFromUCSC(infile, outfile):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses = None
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC(dbhandle, repclasses, outfile)


@follows(mkdir('ucsc.dir'))
@transform(os.path.join(PARAMS["ucsc_dir"],
                        "gbdb",
                        PARAMS["ucsc_database"],
                        "bbi",
                        "*rgMapability*.bw"),
           regex(".*rgMapabilityAlign(\d+)mer.bw"),
           add_inputs(os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"] + ".fasta")),
           r"ucsc.dir/mapability_\1.bed.gz")
def buildMapableRegions(infiles, outfile):
    '''build bed files with mappable regions.

    Convert bigwig tracks with mappability information to a
    bed-formatted file that contains only mappable regions of the
    genome.

    A mapable region is more permissive than a mapable position.

    This method assumes that files use the ``CRG Alignability
    tracks``.

    UCSC says:

      The CRG Alignability tracks display how uniquely k-mer sequences
      align to a region of the genome. To generate the data, the
      GEM-mappability program has been employed. The method is
      equivalent to mapping sliding windows of k-mers (where k has been
      set to 36, 40, 50, 75 or 100 nts to produce these tracks) back to
      the genome using the GEM mapper aligner (up to 2 mismatches were
      allowed in this case). For each window, a mapability score was
      computed (S = 1/(number of matches found in the genome): S=1 means
      one match in the genome, S=0.5 is two matches in the genome, and
      so on). The CRG Alignability tracks were generated independently
      of the ENCODE project, in the framework of the GEM (GEnome
      Multitool) project.

    For the purpose of these tracks, a region is defined to be un-mapable
    if its maximum mapability score is less than 0.5.

    Unmapable regions that are less than half kmer size are mapable, as
    reads from the left/right mapable positions will extend into the region.
    '''

    infile, fastafile = infiles
    fasta = IndexedFasta.IndexedFasta(P.snip(fastafile, ".fasta"))
    contigs = fasta.getContigSizes(with_synonyms=False)

    kmersize = int(re.search(".*Align(\d+)mer.bw", infile).groups()[0])

    E.info("creating mapable regions bed files for kmer size of %i" % kmersize)

    max_distance = kmersize // 2

    f = open(infile)
    bw = BigWigFile(file=f)

    def _iter_mapable_regions(bw, contig, size):

        min_score = PARAMS["ucsc_min_mappability"]

        # there is no iterator access, results are returned as list
        # thus proceed window-wise in 10Mb windows
        window_size = 10000000
        last_start, start = None, None

        for window_start in xrange(0, size, window_size):
            values = bw.get(contig, window_start, window_start + window_size)
            if values is None:
                continue

            for this_start, this_end, value in values:
                if value < min_score:
                    if start:
                        yield start, this_start
                    start = None
                else:
                    if start is None:
                        start = this_start

        if start is not None:
            yield start, this_end

    outf = IOTools.openFile(outfile, "w")

    for contig, size in contigs.iteritems():

        last_start, last_end = None, None
        for start, end in _iter_mapable_regions(bw, contig, size):
            if last_start is None:
                last_start, last_end = start, end
            if start - last_end >= max_distance:
                outf.write("%s\t%i\t%i\n" % (contig, last_start, last_end))
                last_start = start

            last_end = end

        if last_start is not None:
            outf.write("%s\t%i\t%i\n" % (contig, last_start, last_end))

    outf.close()


@transform(buildMapableRegions, suffix(".bed.gz"),
           ".filtered.bed.gz")
def filterMapableRegions(infile, outfile):
    '''remove small windows from a mapability track.

    Too many fragmented regions will cause gat to fail as it
    fragments the workspace into too many individual segments.

    The filtering works by merging all segments that are
    within mapability_merge_distance and removing all those
    that are larger than mapabpility_min_segment_size
    '''

    statement = '''
    mergeBed -i %(infile)s -d %(mapability_merge_distance)i
    | awk '$3 - $2 >= %(mapability_min_segment_size)i'
    | gzip
    > %(outfile)s
    '''

    P.run()


# ---------------------------------------------------------------
# GWAS data
if PARAMS["genome"].startswith("hg"):

    @follows(mkdir('gwas.dir'))
    @merge(None, "gwas.dir/gwascatalog.txt")
    def downloadGWASCatalog(infile, outfile):
        '''download GWAS catalog.'''

        if os.path.exists(outfile):
            os.path.remove(outfile)
        statement = '''wget http://www.genome.gov/admin/gwascatalog.txt
        -O %(outfile)s'''
        P.run()

    @merge(downloadGWASCatalog, PARAMS["interface_gwas_catalog_bed"])
    def buildGWASCatalogTracks(infile, outfile):

        reader = csv.DictReader(IOTools.openFile(infile), dialect="excel-tab")

        tracks = collections.defaultdict(lambda: collections.defaultdict(list))

        fasta = IndexedFasta.IndexedFasta(
            os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"))
        contigsizes = fasta.getContigSizes()
        c = E.Counter()

        for row in reader:
            c.input += 1
            contig, pos, snp, disease = row['Chr_id'], row[
                'Chr_pos'], row['SNPs'], row['Disease/Trait']

            # skip SNPs on undefined contigs
            if contig not in contigsizes:
                c.no_contig += 1
                continue

            if snp == "NR":
                c.skipped += 1
                continue

            if pos == "":
                c.no_pos += 1
                continue

            # translate chr23 to X
            if contig == "23":
                contig = "X"

            contig = "chr%s" % contig

            try:
                tracks[disease][contig].append(int(pos))
            except ValueError:
                print row
            c.output += 1

        E.info(c)

        extension = PARAMS["gwas_extension"]

        c = E.Counter()
        outf = IOTools.openFile(outfile, "w")
        for disease, pp in tracks.iteritems():

            for contig, positions in pp.iteritems():
                contigsize = contigsizes[contig]
                regions = [(max(0, x - extension),
                            min(contigsize, x + extension))
                           for x in positions]

                regions = Intervals.combine(regions)
                c[disease] += len(regions)

                for start, end in regions:
                    outf.write("%s\t%i\t%i\t%s\n" %
                               (contig, start, end, disease))

        outf.close()

        outf = IOTools.openFile(outfile + ".log", "w")
        outf.write("category\tcounts\n%s\n" % c.asTable())
        outf.close()

    @follows(mkdir('gwas.dir'))
    @merge(None, "gwas.dir/gwas_distild.log")
    def downloadDistiLD(infile, outfile):
        '''download GWAS data from distild database.'''

        track = P.snip(outfile, ".log")
        of = track + "_snps.tsv.gz"
        if os.path.exists(of):
            os.path.remove(of)
        statement = \
            '''wget http://distild.jensenlab.org/snps.tsv.gz
            -O %(of)s'''
        P.run()

        of = track + "_lds.tsv.gz"
        if os.path.exists(of):
            os.path.remove(of)
        statement = \
            '''wget http://distild.jensenlab.org/lds.tsv.gz
            -O %(of)s'''
        P.run()

        P.touch(outfile)

    @merge(downloadDistiLD, PARAMS["interface_gwas_distild_bed"])
    def buildDistiLDTracks(infile, outfile):
        '''build bed tracks from DistiLD database.'''

        track = P.snip(infile, ".log")
        intervals = []

        fasta = IndexedFasta.IndexedFasta(
            os.path.join(PARAMS["genome_dir"],
                         PARAMS["genome"] + ".fasta"))
        contigsizes = fasta.getContigSizes()

        c = E.Counter()
        for line in IOTools.openFile(track + "_snps.tsv.gz"):
            pubmed_id, rs, pvalue, block, ensgenes, short, icd10 = line[
                :-1].split("\t")
            c.input += 1
            try:
                contig, start, end = re.match(
                    "(\S+):(\d+)-(\d+)", block).groups()
            except AttributeError:
                E.warn("parsing error for %s" % block)
                c.errors += 1
                continue

            # skip SNPs on undefined contigs
            if contig not in contigsizes:
                c.no_contig += 1
                continue

            intervals.append((contig, int(start), int(end), short))
            c.parsed += 1

        intervals.sort()
        outf = IOTools.openFile(outfile, "w")
        cc = E.Counter()
        for k, x in itertools.groupby(intervals, key=lambda x: x):
            outf.write("%s\t%i\t%i\t%s\n" % k)
            c.output += 1
            cc[k[3]] += 1
        outf.close()
        E.info(c)

        outf = IOTools.openFile(outfile + ".log", "w")
        outf.write("category\tcounts\n%s\n" % cc.asTable())
        outf.close()

    @follows(buildGWASCatalogTracks, buildDistiLDTracks)
    def _gwas():
        pass
else:
    @files(((None, None),))
    def _gwas():
        pass


# ---------------------------------------------------------------
# Ontologies
@follows(mkdir('ontologies.dir'))
@files([(None, PARAMS["interface_go_ensembl"]), ])
def createGO(infile, outfile):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGOFromENSEMBL(infile, outfile)


@transform(createGO,
           regex("(.*)"),
           PARAMS["interface_goslim_ensembl"])
def createGOSlim(infile, outfile):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGOSlimFromENSEMBL(infile, outfile)


@transform((createGO, createGOSlim),
           suffix(".tsv.gz"),
           r"\1_assignments.load")
def loadGOAssignments(infile, outfile):

    table = P.toTable(outfile)
    statement = '''
    zcat < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --table=%(table)s
              --index=gene_id
              --index=go_id
    > %(outfile)s
    '''
    P.run()


@transform(createGO, suffix(".tsv.gz"), ".paths")
def buildGOPaths(infile, outfile):
    '''compute a file with paths of each GO term to the ancestral node.'''
    infile = P.snip(infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOPaths(infile, outfile)


@transform(createGO, suffix(".tsv.gz"), ".desc.tsv")
def buildGOTable(infile, outfile):
    '''build a simple table with GO descriptions in obo.'''
    infile = P.snip(infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOTable(infile, outfile)


@transform(buildGOTable, suffix(".tsv"), ".load")
def loadGOTable(infile, outfile):
    '''load GO descriptions into database.'''
    P.load(infile, outfile)


@follows(mkdir('ontologies.dir'),
         downloadTranscriptInformation, loadGOAssignments)
@files([(None, PARAMS["interface_go_geneontology"]), ])
def createGOFromGeneOntology(infile, outfile):
    '''build GO assignments from GeneOntology.org'''
    PipelineGO.createGOFromGeneOntology(infile, outfile)


@transform(createGOFromGeneOntology,
           suffix(".tsv.gz"),
           add_inputs(buildGOPaths),
           PARAMS["interface_go_geneontology_imputed"])
def imputeGO(infiles, outfile):
    '''imput ancestral GO terms for each gene based on
    derived GO terms.
    '''
    PipelineGO.imputeGO(infiles[0], infiles[1], outfile)


@follows(mkdir('ontologies.dir'))
@files(None, PARAMS['interface_kegg'])
def importKEGGAssignments(infile, outfile):
    ''' import the KEGG annotations from the R KEGG.db
    annotations package. Note that since KEGG is no longer
    publically availible, this is not up-to-date and maybe removed
    from bioconductor in future releases '''

    biomart_dataset = PARAMS["KEGG_dataset"]
    mart = PARAMS["KEGG_mart"]
    host = PARAMS["KEGG_host"]

    PipelineKEGG.importKEGGAssignments(outfile, mart, host, biomart_dataset)


@transform(importKEGGAssignments, suffix(".tsv.gz"), "_assignments.load")
def loadKEGGAssignments(infile, outfile):

    P.load(infile, outfile, options="-i gene_id -i kegg_id")


# ---------------------------------------------------------------
# Enrichment analysis
@follows(mkdir('enrichment.dir'))
@files(buildGeneSet, PARAMS['interface_annotation_gff'])
def annotateGenome(infile, outfile):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes.
    Processed_transcripts tend to cover larger genomic regions
    and often overlap between adjacent protein coding genes.

    In case of overlapping genes, only take the longest
    (in genomic coordinates).

    Genes not on UCSC contigs are removed.
    '''
    PipelineGeneset.annotateGenome(infile,
                                   outfile,
                                   only_proteincoding=True)


@follows(mkdir('enrichment.dir'))
@files(buildGeneSet, PARAMS['interface_genestructure_gff'])
def annotateGeneStructure(infile, outfile):
    '''annotate genome with gene structures.

    Only considers protein coding genes.

    Processed_transcripts tend to cover larger genomic regions
    and often overlap between adjacent protein coding genes.

    In case of overlapping genes, only take the longest
    (in genomic coordinates).

    Genes not on UCSC contigs are removed.
    '''
    PipelineGeneset.annotateGeneStructure(infile,
                                          outfile,
                                          only_proteincoding=True)


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_territories_gff"])
def buildGeneTerritories(infile, outfile):
    '''build gene territories from protein coding genes.'''

    statement = '''
    gunzip < %(infile)s
    | awk '$2 == "protein_coding"'
    | python %(scriptsdir)s/gtf2gtf.py
    --sort=gene
    | python %(scriptsdir)s/gtf2gtf.py
    --merge-transcripts --with-utr
    | python %(scriptsdir)s/gtf2gtf.py
    --sort=position
    | python %(scriptsdir)s/gtf2gff.py
          --genome-file=%(genome_dir)s/%(genome)s
          --log=%(outfile)s.log
          --radius=%(enrichment_territories_radius)s
          --method=territories
    | python %(scriptsdir)s/gtf2gtf.py
    --filter=longest-gene --log=%(outfile)s.log
    | gzip
    > %(outfile)s '''

    P.run()


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_tssterritories_gff"])
def buildTSSTerritories(infile, outfile):
    '''build gene territories from protein coding genes.'''

    statement = '''
    gunzip < %(infile)s
    | awk '$2 == "protein_coding"'
    | python %(scriptsdir)s/gtf2gtf.py
    --filter=representative-transcript --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py --sort=position
    | python %(scriptsdir)s/gtf2gff.py
          --genome-file=%(genome_dir)s/%(genome)s
          --log=%(outfile)s.log
          --radius=%(enrichment_territories_radius)s
          --method=tss-territories
    | python %(scriptsdir)s/gtf2gtf.py
    --sort=gene+transcript --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
    --filter=longest-gene --log=%(outfile)s.log
    | gzip
    > %(outfile)s '''

    P.run()


@follows(mkdir('enrichment.dir'))
@merge(buildFlatGeneSet, PARAMS["interface_greatdomains_gff"])
def buildGREATRegulatoryDomains(infile, outfile):
    '''build gene territories from protein coding genes.'''

    statement = '''
    gunzip < %(infile)s
    | awk '$2 == "protein_coding"'
    | python %(scriptsdir)s/gtf2gtf.py
    --filter=representative-transcript
    --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gff.py
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log
    --method=great-domains
    --radius=%(enrichment_great_radius)s
    --upstream=%(enrichment_great_upstream)i
    --downstream=%(enrichment_great_downstream)i
    | gzip
    > %(outfile)s '''

    P.run()


@follows(mkdir('enrichment.dir'))
@merge((importRepeatsFromUCSC,
        importRNAAnnotationFromUCSC,
        buildGeneSet,
        buildFlatGeneSet,
        importCpGIslandsFromUCSC,
        createGO,
        ),
       PARAMS["interface_genomic_context_bed"])
def buildGenomicContext(infiles, outfile):
    '''build a file with genomic context.

    The output is a bed formatted file, annotating genomic segments
    according to whether they are any of the ENSEMBL annotations.

    It also adds the RNA and repeats annotations from the UCSC.

    The annotations can be partially or fully overlapping.

    Adjacent features (less than 10 bp apart) of the same type are merged.
    '''
    PipelineGeneset.buildGenomicContext(infiles, outfile)


@transform(buildGenomicContext, suffix(".bed.gz"), ".tsv")
def buildGenomicContextStats(infile, outfile):
    '''analysis overlap of genomic contexts.'''

    tmpdir = P.getTempDir(".")

    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/split_file.py
        --pattern-output=%(tmpdir)s/%%s.bed
        --column=4
    > %(outfile)s.log
    '''

    P.run()

    statement = '''
    python %(scriptsdir)s/diff_bed.py
       %(tmpdir)s/*.bed
    > %(outfile)s
    '''
    P.run()

    shutil.rmtree(tmpdir)


@follows(mkdir('enrichment.dir'))
@merge((buildGeneTerritories, loadGOAssignments),
       (PARAMS["interface_genomic_function_bed"],
        PARAMS["interface_genomic_function_tsv"],
        ))
def buildGenomicFunctionalAnnotation(infiles, outfiles):
    '''output a bed file with genomic regions with functional annotations.

    Each bed entry is a gene territory. Bed entries are labeled
    by functional annotations associated with a gene.

    Ambiguities in territories are resolved by outputting
    annotations for all genes within a territory.

    The output file contains annotations for both GO and GOSlim. These
    are prefixed by ``go:`` and ``goslim:``.
    '''

    territories_file = infiles[0]

    dbh = connect()

    PipelineGeneset.buildGenomicFunctionalAnnotation(
        territories_file, dbh, outfiles)


@files((buildGeneSet,
        buildPeptideFasta),
       PARAMS["interface_pseudogenes_gtf"])
def buildPseudogenes(infile, outfile):
    '''build set of pseudogenes.'''
    dbh = connect()
    PipelineGeneset.buildPseudogenes(infile, outfile, dbh)


@follows(mkdir('geneset.dir'))
@files((None,),
       PARAMS["interface_numts_gtf"])
def buildNUMTs(infile, outfile):
    '''build list of NUMTs.'''
    PipelineGeneset.buildNUMTs(infile, outfile)

# --------------------------------------------
# Below is a collection of functions that are
# currently inactivated.


if 0:
    ############################################################
    ############################################################
    ############################################################
    # get UCSC tables
    ############################################################
    def getUCSCTracks(infile=PARAMS["filename_ucsc_encode"]):
        '''return a list of UCSC tracks from infile.'''
        tables = []
        with open(infile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                tablename = line[:-1].strip()
                if tablename == "":
                    continue
                tables.append(tablename)
        return tables

    ############################################################
    ############################################################
    ############################################################
    # import UCSC encode tracks
    ############################################################
    @posttask(touch_file("ucsc_encode.import"))
    @files(PARAMS["filename_ucsc_encode"], "ucsc_encode.import")
    def importUCSCEncodeTracks(infile, outfile):

        statement = '''
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu
                 -A -B -e "SELECT * FROM %(tablename)s" %(ucsc_database)s |\
        python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --table=%(tablename)s \
        >> %(outfile)s

        '''

        dbhandle = sqlite3.connect(PARAMS["database"])

        cc = dbhandle.cursor()
        tables = set(
            [x[0] for x in cc.executewait(
                dbhandle,
                "SELECT name FROM sqlite_master WHERE type='table'")])
        cc.close()

        for tablename in getUCSCTracks(infile):
            if tablename in tables:
                E.info("skipping %(tablename)s - already exists" % locals())
                continue

            E.info("importing %(tablename)s" % locals())
            P.run()

    ############################################################
    ############################################################
    ############################################################
    # export UCSC encode tracks as bed
    ############################################################
    @transform(importUCSCEncodeTracks, suffix(".import"), ".bed")
    def exportUCSCEncodeTracks(infile, outfile):

        dbhandle = sqlite3.connect(PARAMS["database"])

        outs = open(outfile, "w")
        for tablename in getUCSCTracks():
            outs.write("track name=%s\n" % tablename)

            cc = dbhandle.cursor()
            statement = """SELECT chrom, chrostart, chroend FROM %s
            ORDER by chrom, chrostart""" % (
                tablename)
            cc.executewait(dbhandle, statement)
            for contig, start, end in cc:
                outs.write("%s\t%i\t%i\n" % (contig, start, end))
        outs.close()


# ############################################################
# ############################################################
# ############################################################
# # Mappability
# @files(os.path.join(PARAMS["gem_dir"],
#                     PARAMS["genome"] + ".gem"),
#        PARAMS["genome"] + ".mappability")
# def calculateMappability(infile, outfile):
#     '''Calculate mappability using GEM '''
#     index = P.snip(infile, ".gem")
#     job_threads = PARAMS["gem_threads"]
#     statement = '''gem-mappability
#     -t %(gem_threads)s -m %(gem_mismatches)s
#     --max-indel-length %(gem_max_indel_length)s
#     -l %(gem_window_size)s
#     -I %(index)s -o %(outfile)s '''
#     P.run()


# @transform(calculateMappability, suffix(".mappability"),
#            ".mappability.count")
# def countMappableBases(infile, outfile):
#     '''Count mappable bases in genome'''
#     statement = '''cat %(infile)s | tr -cd ! | wc -c > %(outfile)s'''
#     P.run()


# @transform(countMappableBases, suffix(".count"), ".count.load")
# def loadMappableBases(infile, outfile):
#     '''load count of mappable bases in genome'''
#     header = "total_mappable_bases"
#     statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
#                       --table=total_mappable_bases
#                       --header=%(header)s
#                    > %(outfile)s '''
#     P.run()

# ###################################################################


# @transform(calculateMappability, suffix(".mappability"), ".split.log")
# def splitMappabiliyFileByContig(infile, outfile):
#     '''Count mappable bases in genome'''
#     track = P.snip(os.path.basename(infile), ".mappability")
#     statement = '''mkdir contigs;
#     csplit -k -f contigs/contig %(infile)s '/^~[a-zA-Z]/' {100000}
#     > %(outfile)s;
#     rm contigs/contig00;'''
#     P.run()

# ###################################################################


# @follows(splitMappabiliyFileByContig)
# @merge("contigs/contig*", PARAMS["genome"] + "_mappability_per_contig.tsv")
# def countMappableBasesPerContig(infiles, outfile):
#     '''Count mappable bases for each contig'''
#     for infile in infiles:
#         statement = '''grep '~' %(infile)s
#         | sed s/~//g >> %(outfile)s; cat %(infile)s | tr -cd !
#         | wc -c >> %(outfile)s'''
#         P.run()

#     statement = '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s;'''
#     P.run()

# ###################################################################


# @transform(countMappableBasesPerContig, suffix(".tsv"), ".tsv.load")
# def loadMappableBasesPerContig(infile, outfile):
#     '''load count of mappable bases per contig '''
#     header = "contig,mappable_bases"
#     statement = '''cat %(infile)s
#     | python %(scriptsdir)s/csv2db.py
#     --table=mappable_bases_per_contig
#     --header=%(header)s
#     > %(outfile)s '''
#     P.run()

@transform("*/*.gff.gz",
           suffix(".gff.gz"),
           ".gffsummary.tsv.gz")
def buildGFFSummary(infile, outfile):
    '''summarize genomic coverage of gff file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2coverage.py
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.bed.gz",
           suffix(".bed.gz"),
           ".bedsummary.tsv.gz")
def buildBedSummary(infile, outfile):
    '''summarize genomic coverage of bed file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/bed2stats.py
    --per-contig
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/genomic_context.bed.gz",
           suffix(".bed.gz"),
           ".bednamesummary.tsv.gz")
def buildBedNameSummary(infile, outfile):
    '''summarize genomic coverage of bed file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/bed2stats.py
    --per-name
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gtf.gz",
           suffix(".gtf.gz"),
           ".gtfsummary.tsv.gz")
def buildGTFSummary(infile, outfile):
    '''summarize genomic coverage of bed file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2coverage.py
    --genome-file=%(genome_dir)s/%(genome)s
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gtf.gz",
           suffix(".gtf.gz"),
           ".gtfstats.tsv.gz")
def buildGTFStats(infile, outfile):
    '''summarize genomic coverage of bed file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2stats.py
    --is-gtf
    | gzip > %(outfile)s
    '''
    P.run()


@transform("*/*.gff.gz",
           suffix(".gff.gz"),
           ".gffstats.tsv.gz")
def buildGFFStats(infile, outfile):
    '''summarize genomic coverage of bed file.'''
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2stats.py
    | gzip > %(outfile)s
    '''
    P.run()


@merge(buildGTFStats, 'gtf_stats.load')
def loadGTFStats(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.*).tsv.gz")


@merge(buildGFFStats, 'gff_stats.load')
def loadGFFStats(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.*).tsv.gz")


@transform((buildGFFSummary,
            buildBedSummary,
            buildBedNameSummary,
            buildGTFSummary),
           suffix(".tsv.gz"),
           ".load")
def loadIntervalSummary(infile, outfile):
    P.load(infile, outfile, options='--allow-empty')


##################################################################
# Primary targets
@follows(buildContigSizes,
         buildContigBed,
         buildUngappedContigBed,
         loadGenomeInformation,
         buildGenomeGCProfile,
         buildCpGBed)
def assembly():
    pass


@follows(buildGeneSet,
         loadTranscripts,
         downloadTranscriptInformation,
         loadGeneStats,
         loadCDSStats,
         loadExonStats,
         loadGeneInformation,
         downloadEntrezToEnsembl,
         downloadTranscriptSynonyms,
         buildExonTranscripts,
         buildCodingExonTranscripts,
         buildNonCodingExonTranscripts,
         buildPseudogenes,
         buildNUMTs,
         buildSelenoList,
         )
def ensembl():
    pass


@follows(buildPeptideFasta,
         buildCDSFasta,
         buildCDNAFasta)
def fasta():
    pass


@follows(buildTranscriptRegions,
         buildTranscriptTSS,
         buildTranscriptTTS,
         buildGeneRegions,
         buildGeneTSS,
         buildGeneTTS,
         buildGeneTSSInterval,
         buildIntergenicRegions)
def geneset():
    pass


@follows(importRepeatsFromUCSC,
         importRNAAnnotationFromUCSC,
         importCpGIslandsFromUCSC,
         loadRepeats,
         countTotalRepeatLength)
def ucsc():
    pass


@follows(buildGeneTerritories,
         buildTSSTerritories,
         buildGREATRegulatoryDomains,
         annotateGeneStructure,
         annotateGenome,
         buildGenomicContext,
         buildGenomicContextStats,
         buildGenomicFunctionalAnnotation)
def enrichment():
    pass


@follows(loadGOAssignments,
         loadKEGGAssignments)
def ontologies():
    pass


@follows(_gwas)
def gwas():
    pass


@follows(loadGTFStats,
         loadGFFStats,
         loadIntervalSummary)
def summary():
    '''summary targets.'''
    pass


# @follows(calculateMappability, countMappableBases,
#          loadMappableBases, splitMappabiliyFileByContig,
#          countMappableBasesPerContig, loadMappableBasesPerContig)
# def gemMappability():
#     '''Count mappable bases in genome'''
#     pass
# taken out gemMappability as not fully configured


@follows(assembly,
         ensembl,
         ucsc,
         geneset,
         fasta,
         ontologies,
         enrichment,
         gwas)
def full():
    '''build all targets.'''
    pass

###################################################################
###################################################################
###################################################################
# primary targets
###################################################################


@follows(mkdir("report"), summary)
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"), summary)
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
