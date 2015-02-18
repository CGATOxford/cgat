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
gpipe/get_genes.py - get genes from database and print in gff format
===============================================================

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

   python gpipe/get_genes.py --help

Type::

   python gpipe/get_genes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import getopt
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.PredictionParser as PredictionParser
import pgdb
import webbrowser


USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Get genes from database and print in gff format.

Version = $Id: gpipe/get_genes.py 1799 2008-03-28 11:44:19Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-C, --Connect=                  database connection string
-S, --schema=                   schema name
-n, --guess-names               guess table names
-E, --table-exons=              table name with exon information
-P, --table-predictions=        table name with predictions
-R, --table-reference=          table name with reference predictions
-G, --table-genes=              table name with genes
-Q, --table-quality=            table name with quality information
-p, --peptides-fasta-file=                 file with input peptide sequences.
-g, --genome-file=           pattern for filenames with the genomic DNA (FASTA).
-f, --format=                   output format [gff,ucsc,gff-seq]
-l, --from-file                 take set of ids from file
-i, --method=filter --filter-method=                   filter to apply
--filter-region=                filter by region (chr1,1234000,1235000[,offset])
--filter-gene=                  filter by gene
--filter-prediction=            filter by prediction
--add-filter-query              add filter for query_token
--restrict-good-exons           only use good exons
--restrict-selected             only use selected transcripts
-a, --ucsc-assembly=            ucsc: assembly to use


Formats are:

gff: minimal gff format
ucsc: ucsc format for track
gff-seq: extended gff format with sequence


""" % sys.argv[0]

param_long_options = ["verbose=", "help", "boundaries=", "exons=", "peptides=", "format=", "genome-file=",
                      "filter=", "filter-gene=", "filter-region=", "filter-prediction=", "add-filter-query=",
                      "table-predictions=", "table-genes=", "table-reference=", "table-exons=", "table-quality=",
                      "ucsc-assembly=", "schema=", "guess-names", "restrict-good-exons", "from-file",
                      "restrict-selected",
                      "version"]

param_short_options = "v:hb:e:p:g:E:P:G:i:a:S:nl"

# pattern for genomes, %s is substituted for the sbjct_token
param_genome_file = "genome_%s.fasta"
param_format = "gff"
param_loglevel = 0
param_filename_peptides = None
param_filter = ""
param_filter_region = None
param_filter_gene = None
param_filter_query = None
param_filter_prediction = None
param_restrict_good_exons = None

param_tablename_exons = None
param_tablename_predictions = None
param_tablename_genes = None
param_tablename_cds = None
param_tablename_redundant = None

param_schemas = None
param_guess_names = 1

param_connection = "db:andreas"

param_ucsc_dir = "/home/andreas/public_html/ucsc_tracks"
param_ucsc_url_host = "www.genome.ucsc.edu"
param_ucsc_url_path = "/cgi-bin/hgTracks"
param_ucsc_assembly = "canFam1"

param_tablename_reference = None

param_from_file = 0

param_restrict_selected = 0

map_class2color = {'CG': "028,028,100",
                   'SG': "028,028,100",
                   'PG': "092,062,190",
                   'UG': "052,052,190",
                   'PP': "140,038,039",
                   'DP': "190,052,052",
                   'SP': "211,100,100",
                   'UP': "200,255,200",
                   'EP': "200,255,200",
                   'CF': "038,038,140",
                   'UF': "052,052,190",
                   'PF': "100,100,211",
                   'UK': "0,0,0",
                   'default': "0,0,0"}


def ReadContigs(dbhandle, tablename):
    genome_lengths = {}
    statement = "SELECT sbjct_token, size, start FROM %s" % tablename
    cc = dbhandle.cursor()
    try:
        cc.execute(statement)
        result = cc.fetchall()
    except pgdb.DatabaseError, msg:
        print "# query failed with message", msg
        result = None

    if result:
        for x, y, z in result:
            genome_lengths[x] = y, z

    return genome_lengths


def BuildLines(dbhandle, statement, genome_lengths, prefix="", default_color=None):

    c = dbhandle.cursor()
    c.execute(statement)

    if param_loglevel >= 2:
        print "# received %i results." % c.rowcount

    sbjct_token = ""
    sbjct_strand = None
    sbjct_from = 10000000000000000
    sbjct_to = 0

    lines = []

    nmatches = 0

    for line in c.fetchall():

        entry = PredictionParser.PredictionParserEntry()

        entry.FillFromTable(line)

        if not genome_lengths.has_key(entry.mSbjctToken):
            filename_genome = param_genome_file % entry.mSbjctToken
            forward_sequences, reverse_sequences = Genomics.ReadGenomicSequences(
                open(filename_genome, "r"))
            genome_lengths[entry.mSbjctToken] = (
                len(forward_sequences[entry.mSbjctToken]), 0)

        lgenome, offset = genome_lengths[entry.mSbjctToken]

        if param_loglevel >= 4:
            print "# lgenome=%i, offset=%i" % (lgenome, offset)

        # get cds information
        exons = []
        if param_tablename_exons:
            cc = dbhandle.cursor()

            statement = """SELECT exon_from, exon_to, exon_frame, genome_exon_from, genome_exon_to
            FROM %s WHERE prediction_id = %i""" % (param_tablename_exons, entry.mPredictionId,)

            if param_restrict_good_exons:
                statement += " AND is_ok = TRUE"

            try:
                cc.execute(statement)
                result = cc.fetchall()
            except pgdb.DatabaseError, msg:
                print "# query failed with message", msg
                result = []

            exons = result
            cc.close()

        if not exons:
            if entry.mMapPeptide2Genome:
                exons = Genomics.Alignment2ExonBoundaries(entry.mMapPeptide2Genome,
                                                          query_from=entry.mQueryFrom -
                                                          1,
                                                          sbjct_from=entry.mSbjctGenomeFrom,
                                                          add_stop_codon=1)
            else:
                exons = [
                    ("", "", 0, entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo)]

        # select gene id
        if param_tablename_genes:
            cc = dbhandle.cursor()
            statement = """SELECT gene_id
            FROM %s WHERE prediction_id = %i""" % (param_tablename_genes, entry.mPredictionId)

            try:
                cc.execute(statement)
                result = cc.fetchone()
            except pgdb.DatabaseError, msg:
                print "# query failed with message", msg
                result = None

            gene_id = result[0]
            dbhandle.commit()
            cc.close()
        else:
            gene_id = 0

        sbjct_token = entry.mSbjctToken
        sbjct_strand = entry.mSbjctStrand
        if entry.mSbjctStrand == "+":
            sbjct_to = max(entry.mSbjctGenomeTo, sbjct_to)
            sbjct_from = min(entry.mSbjctGenomeFrom, sbjct_from)
        else:
            sbjct_to = max(lgenome - entry.mSbjctGenomeFrom, sbjct_to)
            sbjct_from = min(lgenome - entry.mSbjctGenomeTo, sbjct_from)

        if param_format == "gff":
            if entry.mSbjctStrand == "+":
                for exon in exons:
                    lines.append(string.join(map(str, (entry.mSbjctToken,
                                                       "gipe",
                                                       "exon",
                                                       exon[3], exon[4],
                                                       ".",
                                                       "1",
                                                       exon[2],
                                                       "%s%i_%s" % (prefix, entry.mPredictionId, entry.mQueryToken))),
                                             "\t"))
            else:
                exons.reverse()
                for exon in exons:
                    lines.append(string.join(map(str, (entry.mSbjctToken,
                                                       "CPPipe",
                                                       "exon",
                                                       lgenome -
                                                       exon[
                                                           4], lgenome - exon[3],
                                                       ".",
                                                       "-1",
                                                       exon[2],
                                                       "%s%i_%s" % (prefix, entry.mPredictionId, entry.mQueryToken))),
                                             "\t"))

        elif param_format == "ucsc":

            info_lines = []

            name = "%s_%s-%i-%i-%i" % (param_tablename_predictions,
                                       sbjct_token, sbjct_from, sbjct_to,
                                       entry.mPredictionId)

            filename = "%s/%s" % (param_ucsc_dir, name)

            color = default_color
            quality = "??"

            if param_tablename_quality:
                cc = dbhandle.cursor()
                statement = """SELECT prediction_id,
                is_best_prediction, is_conserved, is_partially_conserved,
                has_frameshift, has_stopcodon, class
                FROM %s WHERE prediction_id = %i""" % (param_tablename_quality, entry.mPredictionId)

                try:
                    cc.execute(statement)
                    result = cc.fetchone()
                except pgdb.DatabaseError, msg:
                    print "# query failed with message", msg
                    result = None

                if result:
                    info_lines.append("## Quality information:")
                    info_lines.append(
                        "is_best_prediction= %s" % str(result[1]))
                    info_lines.append("is_conserved= %s" % str(result[2]))
                    info_lines.append(
                        "is_partially_conserved= %s" % str(result[3]))
                    info_lines.append("has_frameshift= %s" % str(result[4]))
                    info_lines.append("has_stopcodon= %s" % str(result[5]))
                    info_lines.append("class= %s" % str(result[6]))
                    if result[6]:
                        if result[6] in map_class2color:
                            color = map_class2color[result[6]]
                        else:
                            color = map_class2color['default']
                        quality = result[6]
                else:
                    info_lines.append("## no quality information found.")
                dbhandle.commit()
                cc.close()

            if param_tablename_exons:
                cc = dbhandle.cursor()

                statement = """SELECT
                exon_from, exon_to, exon_frame,
                reference_from, reference_to, reference_frame,
                pidentity, psimilarity,
                nframeshifts, ngaps, nstopcodons, is_ok,
                genome_exon_from, genome_exon_to
                FROM %s WHERE prediction_id = %i""" % (param_tablename_exons, entry.mPredictionId)
                try:
                    cc.execute(statement)
                    result = cc.fetchall()
                except pgdb.DatabaseError, msg:
                    print "# query failed with message", msg
                    result = None

                if result:
                    info_lines.append(string.join(("From", "To", "Phase",
                                                   "From", "To", "Phase",
                                                   "Pide", "Psim",
                                                   "NFrame", "NGaps", "NStop", "IsOk",
                                                   "From", "To"), "\t"))
                    for r in result:
                        r[6] = "%5.2f" % r[6]
                        r[7] = "%5.2f" % r[7]
                        info_lines.append(string.join(map(str, r), "\t"))
                else:
                    info_lines.append("## no exon information found.")
                cc.close()

            if info_lines:

                if gene_id:
                    info_lines.append("Gene: %i" % gene_id)
                else:
                    info_lines.append("## no gene found.")

                outfile = open(filename, "w")
                outfile.write(string.join(info_lines, "\n") + "\n")
                outfile.close()

                url = "http://fgu200.anat.ox.ac.uk:8080/ucsc_tracks/%s" % (
                    name)

                lines.append("track name=%s%i_%s description='%s:%s%i_(%i)_%s: pide=%i cov=%i' color=%s visibility=2 url=%s" %
                             (prefix, entry.mPredictionId, entry.mQueryToken,
                              quality, prefix, entry.mPredictionId, gene_id, entry.mQueryToken,
                              entry.mPercentIdentity, entry.mQueryCoverage, color, url))

            else:
                lines.append("track name=%s%i_%s description='%s:%s%i_%s: pide=%i cov=%i' color=%s visibility=2" %
                             (prefix, entry.mPredictionId, entry.mQueryToken,
                              quality, prefix, entry.mPredictionId, entry.mQueryToken,
                              entry.mPercentIdentity, entry.mQueryCoverage, color))

            if entry.mSbjctStrand == "+":
                for exon in exons:
                    lines.append(string.join(map(str, (entry.mSbjctToken,
                                                       "CPPipe",
                                                       "exon",
                                                       exon[
                                                           3] + offset, exon[4] + offset,
                                                       ".",
                                                       entry.mSbjctStrand,
                                                       exon[2],
                                                       "%s%i_%s" % (prefix, entry.mPredictionId, entry.mQueryToken))),
                                             "\t"))
            else:
                exons.reverse()
                for exon in exons:
                    lines.append(string.join(map(str, (entry.mSbjctToken,
                                                       "CPPipe",
                                                       "exon",
                                                       lgenome -
                                                       exon[
                                                           4] + offset, lgenome - exon[3] + offset,
                                                       ".",
                                                       entry.mSbjctStrand,
                                                       exon[2],
                                                       "%s%i_%s" % (prefix, entry.mPredictionId, entry.mQueryToken))),
                                             "\t"))
        elif param_format == "gff-seq":
            if param_tablename_cds:
                cc = dbhandle.cursor()

                statement = """SELECT
                cds_id, genome_from, genome_to, sequence 
                FROM %s WHERE prediction_id = %i""" % (param_tablename_cds, entry.mPredictionId)
                try:
                    cc.execute(statement)
                    result = cc.fetchall()
                except pgdb.DatabaseError, msg:
                    print "# query failed with message", msg
                    result = None
            if result:
                if entry.mSbjctStrand == "+":
                    for exon in result:
                        lines.append(string.join(map(str, (entry.mSbjctToken,
                                                           "gipe",
                                                           "exon",
                                                           exon[2], exon[1],
                                                           "100",
                                                           "1",
                                                           ".",
                                                           exon[0],
                                                           gene_id,
                                                           entry.mPredictionId, exon[3])),
                                                 "\t"))
                else:
                    result.reverse()
                    for exon in result:
                        lines.append(string.join(map(str, (entry.mSbjctToken,
                                                           "gpipe",
                                                           "exon",
                                                           lgenome -
                                                           exon[
                                                               2] + offset, lgenome - exon[1] + offset,
                                                           "100",
                                                           "-1",
                                                           ".",
                                                           exon[0],
                                                           gene_id,
                                                           entry.mPredictionId, exon[3])),
                                                 "\t"))
        else:
            raise "unknown format %s" % param_format

    c.close()

    return lines, sbjct_token, sbjct_strand, sbjct_from, sbjct_to

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-C", "--connection"):
            param_connection = a
        elif o in ("-E", "--table-exons"):
            param_tablename_exons = a
        elif o in ("-P", "--table-predictions"):
            param_tablename_predictions = a
        elif o in ("-R", "--table-reference"):
            param_tablename_reference = a
        elif o in ("-G", "--table-genes"):
            param_tablename_genes = a
        elif o in ("-Q", "--table-quality"):
            param_tablename_quality = a
        elif o in ("-S", "--schema"):
            param_schemas = string.split(a, ",")
        elif o in ("-n", "--guess-names"):
            param_guess_names = 1
        elif o in ("-g", "--genome-file"):
            param_genome_file = a
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-l", "--from-file"):
            param_from_file = 1
        elif o in ("-i", "--method=filter --filter-method"):
            param_filter = a
        elif o == "--filter-region":
            param_filter_region = a
        elif o == "--filter-prediction":
            param_filter_prediction = a
        elif o == "--filter-gene":
            param_filter_gene = a
        elif o == "--add-filter-query":
            param_filter_query = a
        elif o in ("-f", "--format"):
            param_format = a
        elif o in ("-a", "--ucsc-assembly"):
            param_ucsc_assembly = a
        elif o == "--restrict-good-exons":
            param_restrict_good_exons = 1
        elif o == "--restrict-selected":
            param_restrict_selected = 1

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    if param_guess_names and param_schemas[0]:
        schema = param_schemas[0]
        param_tablename_predictions = "%s.overview" % (schema)
        param_tablename_exons = "%s.exons" % (schema)
        param_tablename_quality = "%s.quality" % (schema)
        param_tablename_cds = "%s.cds" % (schema)
        param_tablename_redundant = "%s.redundant" % (schema)
        if not param_tablename_genes:
            param_tablename_genes = "%s.genes" % (schema)

    param_extra_tables = [(param_tablename_predictions, "p")]
    param_extra_conditions = ["TRUE"]

    if param_filter_region:
        data = string.split(param_filter_region, ",")
        if len(data) != 3:
            raise "wrong number of parameters for region: region is chr,from,to"

        param_extra_conditions.append(
            "sbjct_token='%s' AND (OVERLAP(export_sbjct_genome_from, export_sbjct_genome_to, %s, %s) > 0)" %
            (data[0], data[1], data[2]))

        ids = [""]

    elif param_filter_gene:
        param_extra_tables.append((param_tablename_genes, "g"))
        param_extra_conditions.append(
            ("g.gene_id = '%s' AND g.prediction_id = p.prediction_id"))

        if param_from_file:
            infile = open(param_filter_gene, "r")
            ids = []
            for line in infile:
                if line[0] == "#":
                    continue
                ids.append(line[:-1].split("\t")[0])
        else:
            ids = [param_filter_gene]

    elif param_filter_prediction:
        param_extra_conditions.append("p.prediction_id = '%s'")

        if param_from_file:
            infile = open(param_filter_prediction, "r")
            ids = []
            for line in infile:
                if line[0] == "#":
                    continue
                ids.append(line[:-1].split("\t")[0])
        else:
            ids = [param_filter_prediction]

    if param_filter_query:
        param_extra_conditions.append(
            "query_token like '%%%s%%' %s" % (param_filter_query))
        ids = [""]

    if param_restrict_selected:
        param_extra_tables.append((param_tablename_redundant, "r"))
        param_extra_conditions.append(
            "p.prediction_id = r.rep_prediction_id AND r.rep_prediction_id = r.mem_prediction_id")

    print E.GetHeader()
    print E.GetParams()

    dbhandle = pgdb.connect(param_connection)

    # get length of genomic segments
    genome_lengths = ReadContigs(dbhandle, "%s.contigs" % schema)

    lines = []

    nmatches = 1

    # build SQL parts for tables and conditions
    tables = string.join(
        map(lambda x: "%s AS %s" % x, param_extra_tables), ",")
    conditions = string.join(param_extra_conditions, " AND ")

    for id in ids:

        if param_loglevel >= 1:
            print "# processing %i/%i" % (nmatches, len(ids))

        nmatches += 1

        # build SQL statement
        statement = "SELECT p.* FROM %s WHERE %s ORDER BY export_sbjct_genome_from" % (
            tables, conditions)
        if "%s" in statement:
            statement = statement % id

        if param_loglevel >= 3:
            print "# SQL statement=%s" % statement

        new_lines, sbjct_token, sbjct_strand, sbjct_from, sbjct_to = BuildLines(dbhandle,
                                                                                statement,
                                                                                genome_lengths,
                                                                                prefix="q_" +
                                                                                schema +
                                                                                "_",
                                                                                default_color="100,100,100")

        lines += new_lines

        for schema in param_schemas[1:]:

            if sbjct_strand == "+":
                range_from = sbjct_from
                range_to = sbjct_to
            else:
                lgenome = genome_lengths[sbjct_token]
                range_from = lgenome - sbjct_to
                range_to = lgenome - sbjct_from

            if param_guess_names and param_schemas[0]:
                param_tablename_predictions = "%s.overview" % (schema)
                param_tablename_exons = "%s.exons" % (schema)
                param_tablename_genes = "%s.genes" % (schema)
                param_tablename_quality = "%s.quality" % (schema)
                param_tablename_cds = "%s.cds" % (schema)
                statement = """SELECT p.* FROM %s AS p
                WHERE sbjct_token='%s' AND sbjct_strand='%s' AND OVERLAP(sbjct_genome_from, sbjct_genome_to, %s, %s) > 0
                ORDER BY sbjct_genome_from
                """ % (param_tablename_predictions, sbjct_token, sbjct_strand, range_from, range_to )
            elif param_tablename_reference:
                param_tablename_quality = None
                statement = """SELECT p.* FROM %s AS p
                WHERE sbjct_token='%s' AND sbjct_strand='%s' AND OVERLAP(sbjct_genome_from, sbjct_genome_to, %s, %s) > 0
                ORDER BY sbjct_genome_from
                """ % (param_tablename_reference, sbjct_token, sbjct_strand, range_from, range_to )

            ref_lines, ref_sbjct_token, ref_sbjct_strand, ref_sbjct_from, ref_sbjct_to = BuildLines(dbhandle,
                                                                                                    statement,
                                                                                                    genome_lengths,
                                                                                                    prefix=schema +
                                                                                                    "_",
                                                                                                    default_color="0,0,255")
            lines += ref_lines

    if param_format == "ucsc":
        lgenome, offset = genome_lengths[sbjct_token]
        sbjct_from += offset - 100
        sbjct_to += offset + 100

        lines = ["browser position %s:%i-%i" %
                 (sbjct_token, sbjct_from, sbjct_to)] + lines

        name = "%s_%s-%i-%i" % (param_tablename_predictions,
                                sbjct_token, sbjct_from, sbjct_to)
        filename = "%s/%s" % (param_ucsc_dir, name)

        outfile = open(filename, "w")
        outfile.write(string.join(lines, "\n"))
        outfile.close()

        os.chmod(filename, 0666)

        url = "http://%s%s?db=%s&position=%s:%i-%i&hgt.customText=http://fgu203.anat.ox.ac.uk:8080/~andreas/ucsc_tracks/%s" % (
            param_ucsc_url_host, param_ucsc_url_path, param_ucsc_assembly,
            sbjct_token, sbjct_from, sbjct_to, name)

        print url
        print string.join(lines, "\n")

        print "# opening browser window for:"
        print "#", url
        webbrowser.open_new(url)

    else:
        print string.join(lines, "\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
