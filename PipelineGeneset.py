'''utility tasks for dealing with ENSEMBL gene sets.

Most of this tasks take a geneset (.gtf.gz) from ENSEMBL
as input.
'''

import sys, re, os, tempfile, collections, shutil

import Pipeline as P
import Experiment as E

PARAMS = P.getParameters()

############################################################
############################################################
############################################################
def buildGeneRegions( infile, outfile, only_proteincoding = False ):
    '''annotate genomic regions with reference gene set.

    *infile* is an ENSEMBL gtf file.

    In case of overlapping genes, only take the longest (in genomic coordinates).

    Genes not on UCSC contigs are removed.

    Only considers protein coding genes, if ``only_proteincoding`` is set.
    '''

    to_cluster = True

    if only_proteincoding: filter_cmd = ''' awk '$2 == "protein_coding"' '''
    else: filter_cmd = "cat"

    statement = """
            gunzip < %(infile)s |\
            %(filter_cmd)s |\
            %(scriptsdir)s/gff_sort genepos |\
            python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s |\
            python %(scriptsdir)s/gtf2gtf.py --merge-exons --with-utr --log=%(outfile)s.log |\
            python %(scriptsdir)s/gtf2gtf.py --filter=longest-gene --log=%(outfile)s.log |\
            %(scriptsdir)s/gff_sort pos |\
            python %(scriptsdir)s/gtf2gff.py --genome-file=%(genome)s --log=%(outfile)s.log --flank=%(geneset_flank)s > %(outfile)s
        """
    P.run()

############################################################
############################################################
############################################################
def buildProteinCodingGenes( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    *infile* is an ENSEMBL gtf file.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''

    to_cluster = True

    statement = """gunzip 
            < %(infile)s 
            | awk '$2 == "protein_coding"' 
            | %(scriptsdir)s/gff_sort genepos 
            | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s 
            | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log 
            | python %(scriptsdir)s/gtf2gtf.py --filter=longest-gene --log=%(outfile)s.log 
            | awk '$3 == "exon"' 
            | python %(scriptsdir)s/gtf2gtf.py --set-transcript-to-gene --log=%(outfile)s.log 
            | %(scriptsdir)s/gff_sort genepos 
            | gzip
            > %(outfile)s
        """
    P.run()

############################################################
############################################################
############################################################
def importGeneInformation( infile, outfile, only_proteincoding = False ):
    '''import gene information gleaned from the attributes
    in the gene set gtf file.

    *infile* is an ENSEMBL gtf file.

    '''

    table = outfile[:-len(".import")]

    if only_proteincoding: filter_cmd = ''' awk '$2 == "protein_coding"' '''
    else: filter_cmd = "cat"

    statement = '''
    gunzip < %(infile)s |\
    %(filter_cmd)s |\
    %(scriptsdir)s/gff_sort genepos |\
    python %(scriptsdir)s/gtf2tab.py --full --only-attributes -v 0 |\
    python %(toolsdir)s/csv_cut.py --remove exon_number transcript_id transcript_name protein_id |\
    hsort 1 | uniq |\
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --index=gene_name \
              --map=gene_name:str \
              --table=%(table)s \
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
def importTranscriptInformation( infile, outfile,
                                 only_proteincoding = False ):
                                 
    '''import the transcript set.

    *infile* is an ENSEMBL gtf file.
    '''
    to_cluster = True

    table = outfile[:-len(".import")]

    if only_proteincoding: filter_cmd = ''' awk '$2 == "protein_coding"' '''
    else: filter_cmd = "cat"

    statement = '''gunzip 
    < %(infile)s 
    | %(filter_cmd)s 
    | awk '$3 == "CDS"' 
    | %(scriptsdir)s/gff_sort genepos
    | python %(scriptsdir)s/gtf2tab.py --full --only-attributes -v 0
    | python %(toolsdir)s/csv_cut.py --remove exon_number 
    | hsort 1 | uniq 
    | csv2db.py %(csv2db_options)s 
              --index=transcript_id 
              --index=gene_id 
              --index=gene_name 
              --map=transcript_name:str 
              --map=gene_name:str 
              --table=%(table)s 
    > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
def importGeneStats( infile, outfile ):
    '''import gene statistics to database.

    The *infile* is the *outfile* from :meth:`buildGenes`
    '''

    # do not run on cluster - 32/64 bit incompatible.
    # to_cluster = True

    table = outfile[:-len(".import")]

    statement = '''
    gunzip < %(infile)s |\
    python %(scriptsdir)s/gtf2table.py \
          --log=%(outfile)s.log \
          --genome=%(genome)s \
          --counter=position \
          --counter=length \
          --counter=composition-na |\
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --map=gene_id:str \
              --table=%(table)s \
    > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
def buildProteinCodingTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS are used.
    '''

    to_cluster = True

    statement = '''
    gunzip < %(infile)s |\
    awk '$2 == "protein_coding"' |\
    awk '$3 == "CDS"' |\
    python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
    python %(scriptsdir)s/gtf2gtf.py --remove-duplicates=gene --log=%(outfile)s.log > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
def importTranscripts( infile, outfile ):
    '''import the transcript set.'''
    table = outfile[:-len(".import")]
    
    statement = '''
    gunzip < %(infile)s |\
    python %(scriptsdir)s/gtf2tab.py |\
    csv2db.py %(csv2db_options)s \
              --index=transcript_id \
              --index=gene_id \
              --table=%(table)s \
    > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
def importTranscriptStats( infile, outfile ):
    '''import gene statistics to database.

    The *infile* is the *outfile* from :meth:`buildTranscripts`
    '''

    to_cluster = True

    table = outfile[:-len(".import")]

    statement = '''
    gunzip < %(infile)s |\
    python %(scriptsdir)s/gtf2table.py \
          --log=%(outfile)s.log \
          --genome=%(genome)s \
          --reporter=transcripts \
          --counter=position \
          --counter=length \
          --counter=composition-na |\
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --map=gene_id:str \
              --table=%(table)s \
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
def importProteinStats( infile, outfile ):
    '''import protein statistics to database.

    The *infile* is an ENSEMBL peptide file.
    '''

    to_cluster = True

    table = outfile[:-len(".import")]

    statement = '''
    gunzip < %(infile)s |
    python %(scriptsdir)s/analyze_codonbias_shannon.py 
          --log=%(outfile)s
          --type=aa 
          --section=length 
          --section=hid 
          --section=aa 
          --regex-identifier="(\S+)" |
    sed "s/^id/protein_id/" |
    csv2db.py %(csv2db_options)s 
              --index=protein_id 
              --map=protein_id:str 
              --table=%(table)s 
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
def buildPromotorRegions( infile, outfile ):
    '''annotate promotor regions from reference gene set.'''
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=%(promotor_size)s \
                              --genome-file=%(genome)s --log=%(outfile)s.log 
        | gzip 
        > %(outfile)s
    """
    P.run()

############################################################
############################################################
############################################################
def buildTSSRegions( infile, outfile ):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1.
    '''
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome)s --log=%(outfile)s.log > %(outfile)s
    """
    P.run()

############################################################
############################################################
############################################################
def buildOverlapWithEnsembl( infile, outfile, filename_bed ):
    '''compute overlap of genes in ``infile`` with intervals
    in ``filename_bed`` and import into database.

    If ``filename_bed`` has multiple tracks the overlap will
    be computed for each track separately.

    ``infile`` is the output from :meth:`buildGenes`.
    '''

    to_cluster = True
    statement = '''gunzip 
        < %(infile)s 
        | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts 
        | python %(scriptsdir)s/gff2bed.py --is-gtf 
        | python %(scriptsdir)s/bed2graph.py 
            --output=name 
            --log=%(outfile)s.log 
            - %(filename_bed)s 
        > %(outfile)s
    '''
    P.run()
