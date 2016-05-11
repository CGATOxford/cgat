'''
probeset2gene - 
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

From a list of affymetrix probesets, find the
associated genes and transcripts.

A transcript must overlap the mapped location
of a probeset in order to be considered.

Uses the bioconductor mappings.

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import collections

from rpy2.robjects import r as R
import CGAT.Experiment as E


def getProbeset2Gene(database):
    '''build map relating a probeset to an ENSEMBL gene_id'''

    prefix = database[:-len(".db")]
    mapping = prefix + "ENSEMBL"
    R.library(database)

    # map is a Bimap object
    m = R(mapping)

    result = R.toTable(m)

    mapping = collections.defaultdict(list)
    for probeset_id, gene_id in zip(result["probe_id"],
                                    result["ensembl_id"]):
        mapping[probeset_id].append(gene_id)

    E.info("obtained %i mappings: probes=%i, genes=%i" %
           (len(result),
            len(set(result["probe_id"])),
            len(set(result["ensembl_id"]))))
    return mapping


def getProbeset2Location(database="hgu133plus2.db"):
    '''build map with genomic coordinates for each probeset.

    The mapping is not necessarily unique.
    '''

    R.library(database)

    prefix = database[:-len(".db")]

    contigs = dict(R(prefix + "CHRLENGTHS"))

    # map is a Bimap object
    result2start = R.toTable(R(prefix + "CHRLOC"))
    result2end = R.toTable(R(prefix + "CHRLOCEND"))

    mapping = collections.defaultdict(list)

    # make sure order is the same
    assert result2start["probe_id"] == result2end["probe_id"]

    for probeset_id, contig, start, end in zip(result2start["probe_id"],
                                               result2start["Chromosome"],
                                               result2start["start_location"],
                                               result2end["end_location"]):

        if start < 0:
            start = contigs[contig] - start
            end = contigs[contig] - end

        mapping[probeset_id].append((contig, start, end))

    E.info("mappings: probes=%i, contigs=%i" %
           (len(set(result2start["probe_id"])),
            len(set(result2start["Chromosome"])),
            ))

    return mapping


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="bioconductor database to use [default=%default].")

    parser.add_option("-m", "--mapping", dest="database", type="string",
                      help="bioconductor mapping to use [default=%default].")

    parser.add_option(
        "-g", "--gtf-file", dest="filename_gtf", type="string",
        help="filename with the gene set in gtf format [default=%default].")

    parser.set_defaults(
        database="mouse4302.db",
        mapping="ENSEMBL",
        filename_gtf=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    prefix = options.database[:-len(".db")]

    mapping_probeset2gene = prefix + options.mapping
    mapping_probeset2loc = prefix + "CHRLOC"

    probeset2gene = getProbeset2Gene(
        database=options.database,
    )

    probeset2location = getProbeset2Location(
        database=options.database,
    )

    # gtf = GTF.readAndIndex(
    #     GTF.iterator( IOTools.openFile( options.filename_gtf ) ) )

    counts = E.Counter()

    outfile_notfound = open("notfound.table", "w")

    options.stdout.write("probeset_id\tgene_id\tngenes\n")

    for probeset, locations in probeset2location.iteritems():
        counts.probesets += 1
        gene_ids = probeset2gene[probeset]
        if len(gene_ids) == 0:
            counts.notfound += 1
            continue

        for gene_id in gene_ids:
            options.stdout.write(
                "%s\t%s\t%i\n" % (probeset, gene_id, len(gene_ids)))
        counts.output += 1

    E.info("%s" % str(counts))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
