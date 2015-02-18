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
gpipe/assignments2pairs.py - 
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

   python gpipe/assignments2pairs.py --help

Type::

   python gpipe/assignments2pairs.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.PredictionParser as PredictionParser
import CGAT.IndexedFasta as IndexedFasta

USAGE = """python %s [OPTIONS] < assignments > pairs

Version: $Id: gpipe/assignments2pairs.py 2011 2008-07-04 10:40:51Z andreas $

Take assignments of transcripts to regions and
massage boundaries so that it is likely, that any
missing exons will be included in the prediction.

Algorithm:

1. check if there are missing exons. If there are,
add margin to both ends.

Options:
-h, --help                      print this message.
-v, --verbose                   loglevel.
-m, --max-margin=               maximum margin for missing terminal exons
-n, --min-margin=               minimum margin to add in any case.
-d, --default-margin=           margin to use as a default.
-g, --genome-file=              pattern for filenames with the genomic DNA (FASTA).
-p, --peptides-fasta-file=                 file with protein sequences (FASTA).
-s, --suffix=                   suffix for output filenames
-p, --column-prefix=                   prefix for output filenames
-i, --input-format=             input format, valid options are:
                                predictions [default]
                                minimal (query, sbjct_token/strand, range)
                                exons (query, sbjct_token/strand, phase, rank, from/to, from/to)
-f, --format=                   output format, valid options are:
                                fasta: single fasta file
                                single_fasta: individual files for genomic segments.
                                chunks: create chunks of chunk size
--no-sequence                   do not write sequence into files
-c, --chunk=                    chunk size
-u, --clusters=                 filename with clustering information
-k, --offset-key                use offset key as identifier
-t, --conserve-strand           do not switch strand
-o, --is-forward-coordinates       coordinates are forward coordinates
-r, --max-region                maximum predictions per region
-a, --output-filename-pattern=           filename pattern for output. Should contain one %%i
--combine-exons                 combine exons
--filename-previous               filename with previous output (that are to be skipped).
""" % sys.argv[0]

global_outfile = None
global_chunk_size = 0
global_nchunks = 0

global_keys = {}

options = {}

# ------------------------------------------------------------


def WriteEntry(query_token,
               sbjct_token, sbjct_strand,
               peptide_sequence, genomic_fragment,
               margin_sbjct_from=0, margin_sbjct_to=0,
               min_sbjct_from=0, min_sbjct_to=0,
               lgenome=0,
               region_id=None, region_nr=None, region_max_nr=None):

    global global_outfile, global_chunk_size, global_nchunks, global_keys, global_filename, options

    key = "%s_vs_%s_%s_%i_%i" % (query_token,
                                 sbjct_token, sbjct_strand,
                                 margin_sbjct_from, margin_sbjct_to)

    if global_keys.has_key(key):
        print "# eliminated duplicate key %s in region: %s, %s" % (key, str(region_id), str(region_nr))
        return 0

    global_keys[key] = 1

    if options.offset_key:
        export_key = string.join(
            map(str, (sbjct_token, margin_sbjct_from, lgenome - margin_sbjct_to)), "_")
    else:
        export_key = string.join(map(str, (sbjct_token, sbjct_strand,
                                           margin_sbjct_from, margin_sbjct_to,
                                           min_sbjct_from, min_sbjct_to)), " ")

    filename = "%s%s_vs_%s_%s_%i_%i%s" % (options.prefix,
                                          query_token, sbjct_token, sbjct_strand,
                                          margin_sbjct_from, margin_sbjct_to,
                                          options.suffix)

    if options.format == "single_fasta":
        outfile = open(filename, "w")
        outfile.write(">%s\n%s\n" % (export_key, genomic_fragment))
        outfile.close()

    elif options.format == "fasta":
        sys.stdout.write(">%s %s\n%s\n" %
                         (query_token, export_key, genomic_fragment))

    elif options.format == "chunks":
        global_chunk_size += 1

        if global_chunk_size > options.chunk_size:
            # start new file, if either
            # there are no region ids
            # region_nr is at max
            # there has been no outfile
            if (not region_id) or (region_nr == 1) or not global_outfile:
                if global_outfile:
                    global_outfile.close()
                global_filename = options.filename_output_pattern % global_nchunks
                global_nchunks += 1
                global_outfile = open(global_filename, "w")
                global_chunk_size = 1

        if peptide_sequence:
            if region_id:
                global_outfile.write(string.join(map(str, (query_token,
                                                           peptide_sequence,
                                                           sbjct_token,
                                                           sbjct_strand,
                                                           genomic_fragment,
                                                           margin_sbjct_from,
                                                           margin_sbjct_to,
                                                           min_sbjct_from,
                                                           min_sbjct_to,
                                                           region_id, region_nr, region_max_nr)), "\t") + "\n")
            else:
                global_outfile.write(string.join(map(str, (query_token,
                                                           peptide_sequence,
                                                           sbjct_token,
                                                           sbjct_strand,
                                                           genomic_fragment,
                                                           margin_sbjct_from,
                                                           margin_sbjct_to,
                                                           min_sbjct_from,
                                                           min_sbjct_to)), "\t") + "\n")
        else:
            print "# Warning: no peptide sequence for %s" % query_token
            return 0

    elif options.format == "keys":
        print key + "\t" + "\t".join((query_token, sbjct_token, sbjct_strand))

    elif options.format == "dry-run":
        pass

    return 1

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/assignments2pairs.py 2011 2008-07-04 10:40:51Z andreas $", usage=globals()["__doc__"])

    parser.add_option("--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("-s", "--suffix", dest="suffix", type="string",
                      help="")

    parser.add_option("-p", "--column-prefix", dest="prefix", type="string",
                      help="")

    parser.add_option("-a", "--output-filename-pattern", dest="filename_output_pattern", type="string",
                      help="")

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="")

    parser.add_option("-i", "--input-format", dest="input_format", type="string",
                      help="")

    parser.add_option("-u", "--clusters", dest="filename_clusters", type="string",
                      help="")

    parser.add_option("--filename-previous", dest="filename_previous", type="string",
                      help="")

    parser.add_option("-m", "--max-margin", dest="max_margin", type="int",
                      help="")

    parser.add_option("-n", "--min-margin", dest="min_margin", type="int",
                      help="")

    parser.add_option("-d", "--default-margin", dest="default_margin", type="int",
                      help="")

    parser.add_option("-r", "--max-region", dest="max_region_nr", type="int",
                      help="")

    parser.add_option("-c", "--chunk", dest="chunk_size", type="int",
                      help="")

    parser.add_option("-k", "--offset-key", dest="offset_key", action="store_true",
                      help="")

    parser.add_option("-t", "--conserve-strand", dest="conserve_strand", action="store_true",
                      help="")

    parser.add_option("-o", "--is-forward-coordinates", dest="forward_coordinates", action="store_true",
                      help="")

    parser.add_option("--no-sequence", dest="no_sequence", action="store_true",
                      help="")

    parser.add_option("--combine-exons", dest="combine_exons", action="store_true",
                      help="")

    parser.set_defaults(
        # pattern for genomes, %s is substituted for the sbjct_token
        genome_file="genome_%s.fasta",
        filename_peptides=None,
        # margin to add to genomic segments
        max_margin=0,
        min_margin=0,
        default_margin=0,
        offset_key=0,
        chunk_size=100,
        report_step=1000,
        # wheher to combine exons
        combine_exons=False,
        # output format
        format="single_fasta",
        # prefix/suffix for output files
        suffix=".fasta",
        prefix="",
        filename_clusters=None,
        output=None,
        # conserve strand
        conserve_strand=None,
        # input format
        input_format=None,
        forward_coordinates=None,
        # maximum number of predictions per region (0=all)
        max_region_nr=0,
        filename_output_pattern=None,
        # do not write sequences into output
        no_sequence=None,
        # filename with previous results
        filename_previous=None, )

    (options, args) = E.Start(parser)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(1)

    if not options.filename_output_pattern:
        options.filename_output_pattern = options.prefix + \
            "%i" + options.suffix

    # read peptide sequences
    if options.filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(options.filename_peptides, "r"))
    else:
        peptide_sequences = {}

    if options.loglevel >= 1:
        print "# read %i peptide sequences." % len(peptide_sequences)
        sys.stdout.flush()

    # read clustering information
    if options.filename_clusters:
        # Note: if there are no alternative transcripts, map_rep2mem and map_mem2rep will be empty.
        # thus add some dummy variables so that filtering will work.
        map_rep2mem, map_mem2rep = Genomics.ReadMap(
            open(options.filename_clusters, "r"))
        map_rep2mem['dummy'] = ["dummy", ]
        map_mem2rep['dummy'] = "dummy"
    else:
        map_rep2mem, map_mem2rep = {}, {}

    if options.loglevel >= 1:
        print "# read members: mem2rep=%i, rep2mem=%i" % (len(map_mem2rep), len(map_rep2mem))
        sys.stdout.flush()

    map_previous = {}
    # read previous data
    if options.filename_previous:

        entry = PredictionParser.PredictionParserEntry()

        infile = open(options.filename_previous, "r")
        for line in infile:
            if line[0] == "#":
                continue

            if options.input_format == "graph":
                data = line.split("\t")
                (region_id, region_nr, region_max_nr, sbjct_token, sbjct_strand, region_from, region_to,
                 query_token, weight) = data[:9]

                entry.Read("\t".join(data[9:]))

                key = "%s_vs_%s_%s" % (query_token, sbjct_token, sbjct_strand)

                if key not in map_previous:
                    map_previous[key] = [
                        (entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo)]
                else:
                    map_previous[key].append(
                        (entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))

        if options.loglevel >= 1:
            print "# read %i old entries." % (len(map_previous))
            sys.stdout.flush()

    # variables for file numbering
    global_nchunks = 0
    global_chunk_size = options.chunk_size
    global_outfile = None

    # counters of pairs/regions
    npairs = 0
    nregions = 0
    nskipped = 0

    region_id = None
    region_nr = None
    region_max_nr = None
    last_region_id = None

    last_margin_sbjct_from, last_margin_sbjct_to = None, None

    segments = []
    map_query2segments = {}

    entry = PredictionParser.PredictionParserEntry()

    for line in sys.stdin:

        if line[0] == "#":
            continue

        try:

            if options.input_format == "minimal":
                (entry.mQueryToken, entry.mSbjctToken, entry.mSbjctStrand,
                 entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo) = line[:-1].split("\t")
                entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo = int(
                    entry.mSbjctGenomeFrom), int(entry.mSbjctGenomeTo)
            elif options.input_format == "ensembl":
                (dummy,
                 entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                 entry.mSbjctStrand, entry.mSbjctToken,
                 entry.mQueryToken) = line[:-1].split("\t")
                if entry.mSbjctStrand == "1":
                    entry.mSbjctStrand = "+"
                else:
                    entry.mSbjctStrand = "-"
                entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo = int(
                    entry.mSbjctGenomeFrom), int(entry.mSbjctGenomeTo)
            elif options.input_format == "graph":
                data = line.split("\t")
                (region_id, region_nr, region_max_nr, sbjct_token, sbjct_strand, region_from, region_to,
                 query_token, weight) = data[:9]

                entry.Read("\t".join(data[9:]))

                if map_previous:
                    key = "%s_vs_%s_%s" % (query_token,
                                           sbjct_token, sbjct_strand)

                    if key in map_previous:
                        found = False
                        # check for overlap
                        for a, b in map_previous[key]:
                            if min(b, entry.mSbjctGenomeTo) - max(entry.mSbjctGenomeFrom, a) > 0:
                                found = True
                                break
                        if found:
                            nskipped += 1
                            continue

                region_nr, region_max_nr = map(int, (region_nr, region_max_nr))

                if last_region_id != region_id:
                    nregions += 1
                    last_region_id = region_id

                if options.max_region_nr:
                    region_max_nr = min(region_max_nr, options.max_region_nr)
                    if region_nr > options.max_region_nr:
                        continue
            elif options.input_format == "exons":
                (entry.mQueryToken, entry.mSbjctToken, entry.mSbjctStrand, phase, entry.mRank,
                 peptide_from, peptide_to, entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo) = line[:-1].split("\t")
                entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo, entry.mRank = map(
                    int, (entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo, entry.mRank))
            else:
                entry.Read(line)

            if entry.mSbjctStrand == "1":
                entry.mSbjctStrand = "+"
            if entry.mSbjctStrand == "-1":
                entry.mSbjctStrand = "-"

        except ValueError, IndexError:
            print "# Parsing error line: %s" % line[:-1]
            continue

        # increase margin with minimal range
        if options.min_margin:
            min_sbjct_from = max(
                0, entry.mSbjctGenomeFrom - options.min_margin)
            min_sbjct_to = entry.mSbjctGenomeTo + options.min_margin
        else:
            min_sbjct_from = entry.mSbjctGenomeFrom
            min_sbjct_to = entry.mSbjctGenomeTo

        margin_sbjct_from = min_sbjct_from
        margin_sbjct_to = min_sbjct_to

        # increase margin around putative gene region
        if options.default_margin >= 0:
            margin_sbjct_from = max(0, min_sbjct_from - options.default_margin)
            margin_sbjct_to = min_sbjct_to + options.default_margin
        else:
            if entry.mQueryFrom > 0:
                margin_sbjct_from = max(0, min_sbjct_from - options.max_margin)

            if entry.mQueryTo < entry.mQueryLength:
                margin_sbjct_to = min_sbjct_to + options.max_margin

        segments.append([region_id, region_nr, region_max_nr,
                         min_sbjct_from, min_sbjct_to,
                         margin_sbjct_from, margin_sbjct_to,
                         entry.mQueryToken, entry.mSbjctToken, entry.mSbjctStrand])

        if entry.mQueryToken not in map_query2segments:
            map_query2segments[entry.mQueryToken] = []

        map_query2segments[entry.mQueryToken].append(
            [entry.mSbjctToken, entry.mSbjctStrand, margin_sbjct_from, margin_sbjct_to, len(segments) - 1])

    ninput = len(segments)
    nempty, noverlaps, nadjusted, ngenome_notfound, nhigh_overlaps = 0, 0, 0, 0, 0
    nwith_members, nwithout_members, nsame = 0, 0, 0
    unknown = {}
    nunknown = 0

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    for query_token, matches in map_query2segments.items():

        last_m = matches[0]
        for m in matches[1:]:

            # overlap on same chromosome
            if (last_m[0] == m[0] and last_m[1] == m[1] and
                    min(last_m[3], m[3]) - max(last_m[2], m[2]) > 0):

                overlap = min(last_m[3], m[3]) - max(last_m[2], m[2])
                union = max(last_m[3], m[3]) - min(last_m[2], m[2])

                coverage = float(overlap) / union

                noverlaps += 1

                if coverage >= 0.8:
                    nhigh_overlaps += 1
                else:
                    nadjusted += 1
                    segments[last_m[4]][6] -= overlap / 2
                    segments[m[4]][5] += overlap / 2

            last_m = m

    last_query_token = None
    last_sbjct_token = None
    last_sbjct_strand = None
    fragments = []

    if options.loglevel >= 1:
        print "# processing %i segments." % len(segments)
        sys.stdout.flush()

# segments.sort( lambda x,y: cmp( (x[0], x[7], x[8], x[9], , x[5]),
# (y[0], y[7], y[8], y[9], , y[5]) ) )

    # write segments into files
    for segment in segments:

        (region_id, region_nr, region_max_nr,
         min_sbjct_from, min_sbjct_to,
         margin_sbjct_from, margin_sbjct_to,
         query_token, sbjct_token, sbjct_strand) = segment

        try:
            lgenome = fasta.getLength(sbjct_token)
        except KeyError:
            nunknown += 1
            if sbjct_token not in unknown:
                unknown[sbjct_token] = 0
            unknown[sbjct_token] += 1
            continue

        min_sbjct_to = min(min_sbjct_to, lgenome)
        margin_sbjct_to = min(margin_sbjct_to, lgenome)

        if options.forward_coordinates:
            if Genomics.IsNegativeStrand(sbjct_strand):
                margin_sbjct_from, margin_sbjct_to = lgenome - \
                    margin_sbjct_to, lgenome - margin_sbjct_from
                min_sbjct_from, min_sbjct_to = lgenome - \
                    min_sbjct_to, lgenome - min_sbjct_from

        if options.no_sequence:
            fragment = ""
        else:
            # get genomic sequence
            fragment = fasta.getSequence(sbjct_token, sbjct_strand,
                                         margin_sbjct_from, margin_sbjct_to,
                                         as_array=False)

        if peptide_sequences.has_key(query_token):
            peptide_sequence = peptide_sequences[query_token]
        else:
            peptide_sequence = None

        if map_rep2mem:
            if query_token in map_rep2mem:
                nwith_members += 1
                for member_token in map_rep2mem[query_token]:

                    if peptide_sequences.has_key(member_token):
                        peptide_sequence = peptide_sequences[member_token]
                    else:
                        peptide_sequence = None

                    if peptide_sequence == peptide_sequences[query_token]:
                        # same sequence
                        nsame += 1
                        continue

                    npairs += WriteEntry(member_token,
                                         sbjct_token, sbjct_strand,
                                         peptide_sequence, fragment,
                                         margin_sbjct_from, margin_sbjct_to, min_sbjct_from, min_sbjct_to,
                                         lgenome, region_id, region_nr, region_max_nr)
            else:
                nwithout_members += 1

        elif options.combine_exons:
            # if exons are to be comined, keep entries.
            if query_token != last_query_token:
                if last_query_token:
                    fragments.sort()
                    npairs += WriteEntry(last_query_token,
                                         last_sbjct_token,
                                         last_sbjct_strand,
                                         None,
                                         "".join(map(lambda x: x[1], fragments)))
                fragments = []
                last_query_token = query_token
                last_sbjct_token = sbjct_token
                last_sbjct_strand = sbjct_strand
            # save sequence and coordinate. Sort by coordinate later to ensure that
            # exons are in the right order
            fragments.append((min_sbjct_from, fragment))
        else:
            npairs += WriteEntry(query_token,
                                 sbjct_token, sbjct_strand,
                                 peptide_sequence,
                                 "".join(map(lambda x: x[1], fragments)),
                                 margin_sbjct_from, margin_sbjct_to, min_sbjct_from, min_sbjct_to,
                                 lgenome, region_id, region_nr, region_max_nr)

        if options.loglevel >= 1 and (npairs % options.report_step == 0):
            print "# progress: npairs=%i, nregions=%i" % (npairs, nregions)
            sys.stdout.flush()

    if options.combine_exons:
        fragments.sort()
        npairs += WriteEntry(last_query_token,
                             last_sbjct_token,
                             last_sbjct_strand,
                             None,
                             "".join(map(lambda x: x[1], fragments)))

    if global_outfile:
        global_outfile.close()
        if os.path.getsize(global_filename) == 0:
            os.remove(global_filename)

    print "# read=%i, skipped=%i, overlaps=%i, adjusted=%i, empty=%i, genome=%i, high_overlaps=%i, unknown=%i" %\
          (ninput, nskipped, noverlaps, nadjusted, nempty,
           ngenome_notfound, nhigh_overlaps, nunknown)

    if map_rep2mem:
        print "# members: with=%i, without=%i, pairs=%i, same=%i" % (nwith_members, nwithout_members, npairs + nsame, nsame)

    if nunknown:
        print "# unkown:",
        for k, v in unknown.items():
            print "%s=%i," % (k, v),
        print

    print "# written %i pairs" % npairs
    print "# written %i regions" % nregions

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
