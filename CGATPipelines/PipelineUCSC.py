'''PipelineUCSC.py - utility functions for accessing UCSC data
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides methods for accessing ucsc data.

Usage
-----

This module reads the ``[ucsc]`` section of a configuration file
requiring the following variables to be set:

host
    host to connect to

user
    username to connect with

database
    database to use

Documentation
-------------

Code
----

'''

# for UCSC import
import os
import collections
import MySQLdb
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools

# set from calling module
PARAMS = {}

############################################################
############################################################
############################################################
# UCSC tracks
############################################################


def connectToUCSC():
    dbhandle = MySQLdb.Connect(host=PARAMS["ucsc_host"],
                               user=PARAMS["ucsc_user"])

    cc = dbhandle.cursor()
    cc.execute("USE %s " % PARAMS["ucsc_database"])

    return dbhandle


def getRepeatsFromUCSC(dbhandle, repclasses, outfile):
    '''select repeats from UCSC and write to *outfile* in gff format.

    If *repclasses* is None or an empty list, all repeats will be collected.
    '''

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError("could not find any `rmsk` tables")

    # now collect repeats
    tmpfile = P.getTempFile(".")

    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd,
        '.', strand, '.',
        CONCAT('class \\"', repClass, '\\"; family \\"',
        repFamily, '\\"; repName \\"', repName, '\\";')
        FROM %(table)s"""

        if repclasses:
            repclasses_str = ",".join(
                ["'" + x.strip() + "'" for x in repclasses])
            sql += ''' WHERE repClass in (%(repclasses_str)s) ''' % locals()

        sql = sql % locals()

        E.debug("executing sql statement: %s" % sql)
        cc.execute(sql)
        for data in cc.fetchall():
            tmpfile.write("\t".join(map(str, data)) + "\n")

    tmpfile.close()

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = ['''cat %(tmpfilename)s
    | %(scriptsdir)s/gff_sort pos
    | python %(scriptsdir)s/gff2gff.py
    --method=sanitize=genome
    --skip-missing
    --genome-file=%(genome_dir)s/%(genome)s
    --log=%(outfile)s.log ''']

    if PARAMS["ensembl_remove_contigs"]:
        statement.append(
            ''' --remove-contigs="%(ensembl_remove_contigs)s" ''')

    statement.append('''| gzip > %(outfile)s ''')

    statement = " ".join(statement)

    P.run()

    os.unlink(tmpfilename)


def importRefSeqFromUCSC(infile, outfile, remove_duplicates=True):
    '''import gene set from UCSC database
    based on refseq mappings.

    Outputs a gtf-formatted file a la ENSEMBL.

    Depending on *remove_duplicates*, duplicate mappings are either
    removed or kept.

    Matches to chr_random are ignored (as does ENSEMBL).

    Note that this approach does not work as a gene set, as refseq
    maps are not real gene builds and unalignable parts cause
    differences that are not reconcilable.

    '''

    import MySQLdb
    dbhandle = MySQLdb.Connect(host=PARAMS["ucsc_host"],
                               user=PARAMS["ucsc_user"])

    cc = dbhandle.cursor()
    cc.execute("USE %s " % PARAMS["ucsc_database"])

    duplicates = set()

    if remove_duplicates:
        cc.execute("""SELECT name, COUNT(*) AS c FROM refGene
        WHERE chrom NOT LIKE '%_random'
        GROUP BY name HAVING c > 1""")
        duplicates = set([x[0] for x in cc.fetchall()])
        E.info("removing %i duplicates" % len(duplicates))

    # these are forward strand coordinates
    statement = '''
    SELECT gene.name, link.geneName, link.name, gene.name2, product,
    protAcc, chrom, strand, cdsStart, cdsEnd,
    exonCount, exonStarts, exonEnds, exonFrames
    FROM refGene as gene, refLink as link
    WHERE gene.name = link.mrnaAcc
    AND chrom NOT LIKE '%_random'
    ORDER by chrom, cdsStart
    '''

    outf = IOTools.openFile(outfile, "w")

    cc = dbhandle.cursor()
    cc.execute(statement)

    SQLResult = collections.namedtuple(
        'Result',
        '''transcript_id, gene_id, gene_name, gene_id2, description,
        protein_id, contig, strand, start, end,
        nexons, starts, ends, frames''')

    counts = E.Counter()
    counts.duplicates = len(duplicates)

    for r in map(SQLResult._make, cc.fetchall()):

        if r.transcript_id in duplicates:
            continue

        starts = map(int, r.starts.split(",")[:-1])
        ends = map(int, r.ends.split(",")[:-1])
        frames = map(int, r.frames.split(",")[:-1])

        gtf = GTF.Entry()
        gtf.contig = r.contig
        gtf.source = "protein_coding"
        gtf.strand = r.strand
        gtf.gene_id = r.gene_id
        gtf.transcript_id = r.transcript_id
        gtf.addAttribute("protein_id", r.protein_id)
        gtf.addAttribute("transcript_name", r.transcript_id)
        gtf.addAttribute("gene_name", r.gene_name)

        assert len(starts) == len(ends) == len(frames)

        if gtf.strand == "-":
            starts.reverse()
            ends.reverse()
            frames.reverse()

        counts.transcripts += 1
        i = 0
        for start, end, frame in zip(starts, ends, frames):
            gtf.feature = "exon"
            counts.exons += 1
            i += 1
            gtf.addAttribute("exon_number", i)
            # frame of utr exons is set to -1 in UCSC
            gtf.start, gtf.end, gtf.frame = start, end, "."
            outf.write("%s\n" % str(gtf))

            cds_start, cds_end = max(r.start, start), min(r.end, end)
            if cds_start >= cds_end:
                # UTR exons have no CDS
                # do not expect any in UCSC
                continue
            gtf.feature = "CDS"
            # invert the frame
            frame = (3 - frame % 3) % 3
            gtf.start, gtf.end, gtf.frame = cds_start, cds_end, frame
            outf.write("%s\n" % str(gtf))

    outf.close()

    E.info("%s" % str(counts))


#############################################################
#############################################################
#############################################################
##
#############################################################
def getCpGIslandsFromUCSC(dbhandle, outfile):
    '''get CpG islands from UCSC and save as a bed file.

    The name will be set to the UCSC name.
    '''

    cc = dbhandle.cursor()
    table = "cpgIslandExt"
    sql = """SELECT chrom, chromStart, chromEnd, name
    FROM %(table)s ORDER by chrom, chromStart"""
    sql = sql % locals()

    E.debug("executing sql statement: %s" % sql)
    try:
        cc.execute(sql)
        outfile = IOTools.openFile(outfile, "w")
        for data in cc.fetchall():
            outfile.write("\t".join(map(str, data)) + "\n")
        outfile.close()
    except Exception:
        E.warn("Failed to connect to table %s. %s is empty" % (table, outfile))
        P.touch(outfile)

#############################################################
#############################################################
#############################################################
# Methods for setting up a UCSC Track Hub
#############################################################


def readUCSCFile(infile):
    '''read data within a UCSC formatted file.
    returns a list of key,value items
    '''
    result = []
    for line in infile:
        if line.startswith("#"):
            continue
        if line.strip() == "":
            continue
        data = line[:-1].split()
        result.append((data[0], " ".join(data[1:])))
    return result


def writeUCSCFile(outfile, data):
    '''write a hubfile to outfile from data.'''
    for key, value in data:
        outfile.write(" ".join((key, value)) + "\n")


def readTrackFile(infile):
    '''read a track file.

    Returns a list of tracks, each being a dictionary of keywords.
    '''

    data = readUCSCFile(infile)

    def _yielder(data):
        track, block = None, []
        for key, value in data:
            if key == "track":
                if block:
                    yield track, block
                block = []
                track = value
                continue
            block.append((key, value))
    return list(_yielder(data))


def writeTrackFile(outfile, tracks):
    '''write list of *tracks* to track file *outfile*.

    Returns a list of tracks, each being a dictionary of keywords.
    '''

    for track, trackdata in tracks:
        outfile.write(" ".join(("track", track)) + "\n")
        for key, value in trackdata:
            outfile.write(" ".join((key, value)) + "\n")
        outfile.write("\n")
