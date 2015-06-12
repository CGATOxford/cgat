##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
gpipe/analyze_queries.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gpipe/analyze_queries.py --help

Type::

   python gpipe/analyze_queries.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics

import pgdb

""" program $Id: gpipe/analyze_queries.py 2781 2009-09-10 11:33:14Z andreas $

Analyse query predictions.
"""

###############################################################################


def ClusterPeptidesByHid(peptides):
    """cluster peptide sequences by hid."""
    map_cluster2peptide = {}
    map_peptide2cluster = {}

    # cluster peptides by identity
    # (clumsy sort, use hashes for bigger sets)
    for key, sequence in peptides.items():
        h = Genomics.GetHID(sequence)

        if h not in map_cluster2peptide:
            map_cluster2peptide[h] = []

        map_peptide2cluster[key] = h
        map_cluster2peptide[h].append(key)

    return map_cluster2peptide, map_peptide2cluster

###############################################################################


def ClusterPeptidesByFragment(peptides):
    """cluster peptide sequences by identical subsequences.

    Sort sequences by length. Start with longest sequence
    and remove all shorter sequences that are part of it.
    """
    map_cluster2peptide = {}
    map_peptide2cluster = {}

    # cluster peptides by identity
    # (clumsy sort, use hashes for bigger sets)
    keys = peptides.keys()
    keys.sort(lambda x, y: cmp(len(peptides[x]), len(peptides[y])))
    keys.reverse()
    for x in range(len(keys)):
        kx = keys[x]

        if kx in map_peptide2cluster:
            continue

        map_peptide2cluster[kx] = kx
        cluster = [kx]

        s = peptides[keys[x]]

        for y in range(x + 1, len(keys)):
            ky = keys[y]
            if peptides[ky] in s:
                map_peptide2cluster[ky] = kx
                cluster.append(ky)

        map_cluster2peptide[kx] = cluster

    return map_cluster2peptide, map_peptide2cluster

###############################################################################


def CountFoundGenes(result, map_peptide2cluster=None):
    """count found genes and transcripts.

    Counting is done per transcript and gene separately.
    """

    genes = {}
    found_genes = {}
    found_transcripts = {}
    transcripts = {}

    for rep_token, query_token, npredictions in result:
        genes[rep_token] = 1
        if map_peptide2cluster and query_token in map_peptide2cluster:
            t = map_peptide2cluster[query_token]
        else:
            t = query_token
        transcripts[t] = 1
        if npredictions > 0:
            found_transcripts[t] = 1
            found_genes[rep_token] = 1

    return found_genes.keys(), genes.keys(), found_transcripts.keys(), transcripts.keys()

###############################################################################


def GetQueryInfo(dbhandle, genome, options, subset=None):
    """retrive query information from database."""
    if subset or options.filter_quality:
        # apply filter. This means, values need to be recalculated.

        # select all those queries not in previous set and
        # append to this set.
        statement = """
        SELECT rep_token, query_token
        FROM %s.queries AS q
        """ % (genome)

        cc = dbhandle.cursor()
        cc.execute(statement)
        r = cc.fetchall()
        cc.close()

        data = {}
        for x in r:
            data[x[1]] = [x[0], x[1], 0, 0]

        statement = """
        SELECT p.query_token, rep_prediction_id, mem_prediction_id, c.class
        FROM %s.queries AS q, %s.predictions AS p, %s.redundant AS r, %s.quality AS c
        WHERE q.query_token = p.query_token AND 
        r.mem_prediction_id = p.prediction_id AND
        c.prediction_id = p.prediction_id 
        ORDER BY query_token, rep_prediction_id
        """ % (genome, genome, genome, genome)

        cc = dbhandle.cursor()
        cc.execute(statement)
        r = cc.fetchall()
        cc.close()

        if options.loglevel >= 2:
            print "# retrieved %i lines from queries table in %s" % (len(r), genome)
            sys.stdout.flush()

        if subset:
            r = filter(lambda x: str(x[2]) in subset[genome], r)
        if options.filter_quality:
            r = filter(lambda x: str(x[3]) in options.filter_quality, r)

        if options.loglevel >= 2:
            print "# adding up and filtering: %i lines from queries table in %s" % (len(r), genome)
            sys.stdout.flush()

        last_qq = None
        for qq, rep, mem, c in r:
            if qq != last_qq:
                if last_qq:
                    data[last_qq][2] = all
                    data[last_qq][3] = nr
                all, nr = 0, 0
                last_qq = qq
            all += 1
            if rep == mem:
                nr += 1

        data[last_qq][2] = all
        data[last_qq][3] = nr

        r = data.values()
    else:
        statement = """
        SELECT rep_token, query_token, npredictions, nr_npredictions
        FROM %s.queries
        """ % (genome)

        cc = dbhandle.cursor()
        cc.execute(statement)
        r = cc.fetchall()
        cc.close()

        if options.loglevel >= 2:
            print "# retrieved %i lines from queries table in %s" % (len(r), genome)
            sys.stdout.flush()

    return r

# --------------------------------------------------------------------------


def writeListMissed(outfile,
                    output,
                    genomes,
                    options):
    """write list of missed genes/transcripts."""
    options.stdlog.write("# missed: %i rows\n" % (len(output)))

    outfile.write("identifier\tngenomes\tgenomes\n")
    for transcript, missed_genomes in output.items():
        outfile.write(transcript + "\t" + "%i\t" %
                      len(missed_genomes) + "\t".join(missed_genomes) + "\n")

# --------------------------------------------------------------------------


def writeStatsMissed(outfile,
                     output,
                     genomes,
                     options):
    """write statistics of missed genes/transcripts."""
    counts = [0] * (len(genomes) + 1)

    for entry, missed_genomes in output.items():
        counts[len(missed_genomes)] += 1

    outfile.write("ngenomes\tcounts\tpercent\n")
    total = sum(counts)
    for x in range(1, len(genomes)):
        outfile.write("%i\t%i\t%s\n" % (
            x, counts[x], options.format_percent % (100.0 * float(counts[x]) / total)))

# --------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/analyze_queries.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-g", "--genomes", dest="genomes", type="string",
                      help="genomes to analyse.")
    parser.add_option("-i", "--priority", dest="priority", type="string",
                      help="quality priority.")
    parser.add_option("-s", "--method=sort --sort-order", dest="sort", type="string",
                      help="sort order.")
    parser.add_option("-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="filename with template peptide sequences.")
    parser.add_option("-m", "--methods", dest="methods", type="string",
                      help="methods to apply [missed].")
    parser.add_option("-f", "--method=filter --filter-method", dest="filename_filter", type="string",
                      help="filename with schema|prediction_id|gene to use as filter. The prediction_ids are used for filtering.")
    parser.add_option("-q", "--filter-quality", dest="filter_quality", type="string",
                      help="only consider predictions of given qualities.")
    parser.add_option("--pattern-output", dest="pattern_output", type="string",
                      help="output pattern for multiple file output.")
    parser.add_option("--pattern-stats", dest="pattern_stats", type="string",
                      help="output pattern for multiple statistics output.")
    parser.add_option("--outfile-clusters", dest="outfile_clusters", type="string",
                      help="output filename for clusters.")
    parser.add_option("--infile-clusters", dest="infile_clusters", type="string",
                      help="input filename for clusters.")
    parser.add_option("-n", "--non-redundant", dest="non_redundant", action="store_true",
                      help="use non-redundant set for output.")
    parser.add_option("--clustering-method", dest="clustering_method", type="choice",
                      choices=("fragment", "hid"),
                      help="clustering method to use.")

    parser.set_defaults(
        genomes="",
        priority="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        sort="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        methods="missed",
        peptides=None,
        filename_filter=None,
        separator="|",
        filter_quality=None,
        pattern_output="%s",
        pattern_stats=None,
        clustering_method="fragment",
        outfile_clusters=None,
        infile_clusters=None,
        non_redundant=False,
        format_percent="%5.2f",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.filename_peptides:
        peptides = Genomics.ReadPeptideSequences(
            open(options.filename_peptides, "r"))
    else:
        peptides = {}

    if options.genomes:
        options.genomes = options.genomes.split(",")
    if options.priority:
        options.priority = options.priority.split(",")
    if options.methods:
        options.methods = options.methods.split(",")
    if options.sort:
        options.sort = options.sort.split(",")
    if options.filter_quality:
        options.filter_quality = options.filter_quality.split(",")

    subset = {}
    if options.filename_filter:
        data = map(lambda x: x[:-1].split(options.separator)[:3],
                   filter(lambda x: x[0] != "#", open(options.filename_filter, "r").readlines()))
        for s, p, g in data:
            if s not in subset:
                subset[s] = {}
            subset[s][p] = 1

    if len(options.sort) != len(options.priority):
        raise "different number of classes in sort order and priority order"

    dbhandle = pgdb.connect(options.psql_connection)

    # Cluster peptides
    if options.infile_clusters:
        map_peptide2cluster, map_cluster2peptide = IOTools.ReadMap(
            open(options.infile_clusters, "r"), both_directions=True)
    elif peptides:
        if options.clustering_method == "fragment":
            map_cluster2peptide, map_peptide2cluster = ClusterPeptidesByFragment(
                peptides)
        elif options.clustering_method == "hid":
            map_cluster2peptide, map_peptide2cluster = ClusterPeptidesByHid(
                peptides)
    else:
        map_cluster2peptide = {}
        map_peptide2cluster = {}

    if map_cluster2peptide and options.loglevel >= 1:
        options.stdlog.write("# clustering of peptides: %i cluster for %i peptides\n" % (
            len(map_cluster2peptide), len(map_peptide2cluster)))
        sys.stdout.flush()

    if options.outfile_clusters and not options.infile_clusters:
        options.stdlog.write("# writing clusters to %s\n" %
                             options.outfile_clusters)
        outfile = open(options.outfile_clusters, "w")
        for k, v in map_peptide2cluster.items():
            outfile.write("%s\t%s\n" % (k, v))
        outfile.close()

    for method in options.methods:

        if method == "stats":
            # Count number of missed unique genes/transcripts

            headers = ("species",
                       "genes",
                       "found_genes", "missed_genes",
                       "pfound_genes", "pmissed_genes",
                       "nr_found_genes", "nr_missed_genes",
                       "pnr_found_genes", "pnr_missed_genes",
                       "transcripts",
                       "found_transcripts", "missed_transcripts",
                       "pfound_transcripts", "pmissed_transcripts",
                       "nr_found_transcripts", "nr_missed_transcripts",
                       "pnr_found_transcripts", "pnr_missed_transcripts",
                       )

            options.stdout.write("\t".join(headers) + "\n")

            for genome in options.genomes:

                r = GetQueryInfo(dbhandle, genome, options, subset)

                found_genes, genes, found_transcripts, transcripts =\
                    CountFoundGenes(map(lambda x: (x[0], x[1], x[2]), r),
                                    map_peptide2cluster)

                nrfound_genes, nrgenes, nrfound_transcripts, nrtranscripts =\
                    CountFoundGenes(map(lambda x: (x[0], x[1], x[3]), r),
                                    map_peptide2cluster)

                nfound_genes = len(found_genes)
                nfound_transcripts = len(found_transcripts)
                nnrfound_genes = len(nrfound_genes)
                nnrfound_transcripts = len(nrfound_transcripts)

                ngenes = len(genes)
                ntranscripts = len(transcripts)

                if ngenes == 0 or ntranscripts == 0:
                    continue

                f1 = lambda x: 100 * float(x) / ngenes
                f2 = lambda x: 100 * float(x) / ntranscripts
                options.stdout.write("%s\t%i\t%i\t%i\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f\t%i\t%i\t%i\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f\n" % (
                    genome,
                    ngenes,
                    nfound_genes, ngenes - nfound_genes,
                    f1(nfound_genes), f1(ngenes - nfound_genes),
                    nnrfound_genes, ngenes - nnrfound_genes,
                    f1(nnrfound_genes), f1(ngenes - nnrfound_genes),
                    ntranscripts,
                    nfound_transcripts, ntranscripts - nfound_transcripts,
                    f2(nfound_transcripts), f2(
                        ntranscripts - nfound_transcripts),
                    nnrfound_transcripts, ntranscripts - nnrfound_transcripts,
                    f2(nnrfound_transcripts), f2(ntranscripts - nnrfound_transcripts)))

        elif method == "missed":

            headers = ("species",
                       "genes",
                       "transcripts",
                       "missed_genes",
                       "missed_transcripts",
                       "percent_missed_genes",
                       "percent_missed_transcripts",
                       )

            options.stdout.write("\t".join(headers) + "\n")

            all_missed_genes = {}
            all_missed_transcripts = {}

            for genome in options.genomes:

                r = GetQueryInfo(dbhandle, genome, options, subset)

                if options.non_redundant:
                    found_genes, genes, found_transcripts, transcripts =\
                        CountFoundGenes(map(lambda x: (x[0], x[1], x[3]), r),
                                        map_peptide2cluster)
                else:
                    found_genes, genes, found_transcripts, transcripts =\
                        CountFoundGenes(map(lambda x: (x[0], x[1], x[2]), r),
                                        map_peptide2cluster)

                sg = set(genes)
                missed_genes = sg.difference(found_genes)

                for x in missed_genes:
                    if x not in all_missed_genes:
                        all_missed_genes[x] = []
                    all_missed_genes[x].append(genome)

                sm = set(transcripts)
                missed_transcripts = sm.difference(found_transcripts)

                for x in missed_transcripts:
                    if x not in all_missed_transcripts:
                        all_missed_transcripts[x] = []
                    all_missed_transcripts[x].append(genome)

                options.stdout.write("%s\t%i\t%i\t%i\t%i\t%s\t%s\n" % (genome,
                                                                       len(genes),
                                                                       len(transcripts),
                                                                       len(missed_genes),
                                                                       len(missed_transcripts),
                                                                       options.format_percent % (
                                                                           100.0 * float(len(missed_genes)) / len(genes)),
                                                                       options.format_percent % (100.0 * float(len(missed_transcripts)) / len(transcripts))))

            for section in ("genes", "transcripts"):

                if section == "genes":
                    missed = all_missed_genes
                else:
                    missed = all_missed_transcripts

                writeListMissed(open(options.pattern_output % section, "w"),
                                missed, options.genomes, options)

                if options.pattern_stats:
                    outfile = open(options.pattern_stats % section, "w")
                else:
                    outfile = options.stdout
                    outfile.write("# statistics for %s\n" % section)

                writeStatsMissed(outfile, missed, options.genomes, options)

                if outfile != options.stdout:
                    outfile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
