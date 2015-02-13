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
gpipe/liftover_predictions.py - 
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

   python gpipe/liftover_predictions.py --help

Type::

   python gpipe/liftover_predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import bisect
import CGAT.Genomics as Genomics
import alignlib_lite
import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser
from CGAT.Predictor2 import PredictorExonerate

USAGE = """python %s file1 file2

Version: $Id: gpipe/liftover_predictions.py 2781 2009-09-10 11:33:14Z andreas $

map predictions between two assemblies. This procedure here
only accounts for length differences in gap regions (NNN).
It requires a mapping file like it is output from diff_fasta.py.
        
Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.

""" % sys.argv[0]

parser = E.OptionParser(
    version="%prog version: $Id: gpipe/liftover_predictions.py 2781 2009-09-10 11:33:14Z andreas $")

#####################################################################


def ReadOffsets(infile):
    """read offsets from file.

    Format is:
    sbjct_token first last offsets

    in 0 offsets coordinates.
    """

    map_old2new = {}

    for line in infile:
        if line[0] == "#":
            continue
        chr, res_from, res_to, offset = line[:-1].split("\t")[:4]
        res_from, res_to, offset = map(int, (res_from, res_to, offset))
        if chr not in map_old2new:
            map_old2new[chr] = []
        map_old2new[chr].append((res_from, res_to, offset))

    breakpoints = {}
    endpoints = {}
    offsets = {}

    for chr in map_old2new.keys():
        r = map_old2new[chr]
        r.sort()
        breakpoints[chr] = map(lambda x: x[0], r)
        endpoints[chr] = map(lambda x: x[1], r)
        offsets[chr] = map(lambda x: x[2], r)

    return breakpoints, endpoints, offsets

#####################################################################


def GetOffset(value, breakpoints, endpoints, offsets):
    """return offset for a given coordinate."""

    pos = bisect.bisect(endpoints, value)
    if pos == len(breakpoints):
        return offsets[-1]

    return offsets[pos]

#####################################################################


def MapAlignment(entry, map_a2b):
    """map alignment from sbjct_from on sbjct_strand using
    forware alignment map_a2b (which is shifted by 1).

    Note that this procedure won't detect gaps between states, only
    differences in length within states. If there is a gap between states,
    the whole alignment will be meaningless, but there is no way to
    detect this apart from translating to the amino acid sequence
    and comparing it.
    """

    is_positive = entry.mSbjctStrand == "+"

    if is_positive:
        sbjct_pos = entry.mSbjctGenomeFrom + 1
    else:
        # no -1, as it starts on the residue
        sbjct_pos = map_a2b.getRowTo() - entry.mSbjctGenomeFrom

    last_mapped_pos = map_a2b.mapRowToCol(sbjct_pos)

    if last_mapped_pos == 0:
        raise ValueError, "unmappable starting residue %i" % sbjct_pos

    new_alignment = []

    if is_positive:
        entry.mSbjctGenomeFrom = last_mapped_pos - 1
    else:
        entry.mSbjctGenomeFrom = map_a2b.getColTo() - last_mapped_pos

    total_d = 0
    for state, l_query, l_sbjct in entry.mMapPeptide2Genome[:-1]:

        if is_positive:
            sbjct_pos += l_sbjct
        else:
            sbjct_pos -= l_sbjct

        mapped_pos = map_a2b.mapRowToCol(sbjct_pos)

        if mapped_pos == 0:
            for x in 1, 2:
                if map_a2b.mapRowToCol(sbjct_pos + x):
                    sbjct_pos += x
                    mapped_pos = map_a2b.mapRowToCol(sbjct_pos)
                    break
            else:
                raise ValueError, "unmappable residue %i" % sbjct_pos

        d = abs(mapped_pos - last_mapped_pos)
        total_d += d
        new_alignment.append((state, l_query, d))

        last_mapped_pos = mapped_pos

    state, l_query, l_sbjct = entry.mMapPeptide2Genome[-1]

    # process last state, map only to last residue
    if is_positive:
        sbjct_pos += l_sbjct - 1
    else:
        sbjct_pos -= l_sbjct - 1

    mapped_pos = map_a2b.mapRowToCol(sbjct_pos)

    if mapped_pos == 0:
        raise ValueError, "unmappable residue %i" % sbjct_pos

    d = abs(mapped_pos - last_mapped_pos) + 1
    total_d += d

    new_alignment.append((state, l_query, d))

    entry.mSbjctGenomeTo = entry.mSbjctGenomeFrom + total_d

    entry.mMapPeptide2Genome = new_alignment


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-m", "--map-tsv-file", dest="filename_map", type="string",
                      help="filename with mapping information.")
    parser.add_option("-o", "--pattern-old", dest="pattern_old", type="string",
                      help="pattern for mapping new to old identifiers: extract string from old.")
    parser.add_option("-n", "--pattern-new", dest="pattern_new", type="string",
                      help="pattern for mapping new to old identifiers: put string into new.")
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="genome_file.")
    parser.add_option("-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="filename with peptide sequences.")
    parser.add_option("-f", "--input-format", dest="input_format", type="choice",
                      help="format of mapping file", choices=("alignment", "offsets"))
    parser.add_option("-i", "--write-missed", dest="write_missed", type="string",
                      help="write missed identifiers to separate file.")
    parser.add_option("-a", "--genes-tsv-file", dest="filename_genes", type="string",
                      help="filename with gene information.")
    parser.add_option("--filename-old-peptides", dest="filename_old_peptides", type="string",
                      help="filename with old peptide information.")
    parser.add_option("--no-renumber", dest="renumber", action="store_false",
                      help="do not renumber predictions.")
    parser.add_option("--contig-sizes-old", dest="contig_sizes_old", type="string",
                      help="contig sizes for old data.")
    parser.add_option("--contig-sizes-new", dest="contig_sizes_new", type="string",
                      help="contig sizes for new data.")
    parser.add_option("--skip-errors", dest="skip_errors", action="store_true",
                      help="skip entries with errors.")

    parser.set_defaults(
        filename_map=None,
        pattern_old="(.+)",
        pattern_new="%s",
        genome_file=None,
        filename_peptides=None,
        write_missed=None,
        filename_genes=None,
        filename_old_peptides=None,
        renumber=True,
        input_format="alignment",
        contig_sizes_old=None,
        contig_sizes_new=None,
        skip_errors=None
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    predictor = PredictorExonerate()

    # the different mapping criteria
    map_sbjcts = {}
    breakpoints = {}

    ##########################################################################
    map_transcript2gene = {}
    if options.filename_genes:
        infile = open(options.filename_genes, "r")
        for gene, transcript in map(lambda x: x[:-1].split("\t")[:2], filter(lambda x: x[0] != "#", infile.readlines())):
            map_transcript2gene[transcript] = gene
        infile.close()

    ##########################################################################
    peptides = {}
    if options.filename_peptides:
        peptides = Genomics.ReadPeptideSequences(
            open(options.filename_peptides, "r"))
        options.stdlog.write("# read %i peptide sequences.\n" % len(peptides))

    ##########################################################################
    # read old query sequences and compare against new query sequences
    # this can be used to build a map between old and new queries
    query_map_old2new = {}
    if options.filename_old_peptides:
        old_peptides = Genomics.ReadPeptideSequences(
            open(options.filename_old_peptides, "r"))
        options.stdlog.write(
            "# read %i old peptide sequences.\n" % len(old_peptides))
        query_map_old2new, unmappable, unmapped = Genomics.MapSequences(
            old_peptides, peptides)
        options.stdlog.write(
            "# built map: unmappable=%i unmapped=%i.\n" % (len(unmappable), len(unmapped)))
        if options.loglevel >= 2:
            options.stdlog.write("# unmappable: %s.\n" % ";".join(unmappable))
            options.stdlog.write("# unmapped: %s.\n" % ";".join(unmapped))

    ##########################################################################
    # read old/new contig sizes for mapping positive/negative coordinates
    contig_sizes_old = {}
    contig_sizes_new = {}
    if options.contig_sizes_old:
        contig_sizes_old = Genomics.ReadContigSizes(
            open(options.contig_sizes_old, "r"))
    if options.contig_sizes_new:
        contig_sizes_new = Genomics.ReadContigSizes(
            open(options.contig_sizes_new, "r"))

    ##########################################################################
    if options.filename_map:

        infile = open(options.filename_map)
        if options.input_format == "alignments":
            for line in infile:
                if line[0] == "#":
                    continue

                x, old_token, old_from, old_to, old_ali, new_from, new_to, new_ali = line[
                    :-1].split("\t")

                map_sbjcts[old_token] = (old_from, old_ali, new_from, new_ali)

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# read %i alignments.\n" % len(map_sbjcts))

        elif options.input_format == "offsets":
            # input is a list of segments and their offsets.

            breakpoints, endpoints, offsets = ReadOffsets(infile)
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# read breakpoints for %i chromosomes.\n" % len(breakpoints))

        infile.close()

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # end of input section
    ##########################################################################
    ##########################################################################
    ##########################################################################

    rx = re.compile(options.pattern_old)
    last_sbjct_token = None
    ninput = 0
    nerrors = 0
    nerrors_map = 0
    nerrors_inconsistencies = 0
    nerrors_boundaries = 0
    nerrors_translation = 0
    nerrors_inconsequential = 0
    nerrors_realigned = 0
    nmapped = 0
    nfiltered = 0
    naligned = 0
    noutput = 0
    found_transcripts = {}
    nduplicates = 0
    output = {}

    for line in sys.stdin:
        if line[0] == "#":
            continue

        entry = PredictionParser.PredictionParserEntry()

        entry.Read(line)

        ninput += 1
        is_positive = entry.mSbjctStrand == "+"

        is_error = False

        # check if query token is mappable: using sequence map
        if (query_map_old2new and entry.mQueryToken not in query_map_old2new):
            options.stdlog.write("# skipping prediction %i: obsolete query %s\n" % (
                entry.mPredictionId, entry.mQueryToken))
            nfiltered += 1
            continue
        else:
            # check if query token is mappable: using filter
            if (peptides and entry.mQueryToken not in peptides):
                options.stdlog.write("# skipping prediction %i: obsolete query %s\n" % (
                    entry.mPredictionId, entry.mQueryToken))
                nfiltered += 1
                continue

        new_sbjct_token = options.pattern_new % rx.search(
            entry.mSbjctToken).groups()[0]

        #######################################################################
        # Map via alignments
        if entry.mSbjctToken in map_sbjcts:
            nmapped += 1
            if last_sbjct_token != entry.mSbjctToken:
                old_from, old_ali, new_from, new_ali = map_sbjcts[
                    entry.mSbjctToken]
                map_a2b = alignlib_lite.makeAlignmentVector()
                alignlib_lite.AlignmentFormatExplicit(
                    int(old_from), old_ali,
                    int(new_from), new_ali).copy(map_a2b)

            last_sbjct_token = entry.mSbjctToken

            if options.loglevel >= 3:
                print "#", str(entry)
                print "#", map_sbjcts[entry.mSbjctToken]
                sys.stdout.flush()

            old_f, old_t = entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo

            # convert to forward coordinates:
            if is_positive:
                f, t = old_f, old_t
                first_res, last_res = f + 1, t
            else:
                f, t = map_a2b.getRowTo() - old_f, map_a2b.getRowTo() - old_t
                first_res, last_res = f, t + 1

            # map first and last residues
            mfirst_res = map_a2b.mapRowToCol(first_res)
            mlast_res = map_a2b.mapRowToCol(last_res)

            if (mfirst_res == 0 and old_f != 0) or (mlast_res == 0 and old_t != map_a2b.getRowTo()):

                options.stderr.write("# mapping not possible for prediction %i on %s %s:%i-%i -> %i-%i -> %i-%i -> %i-%i -> %i-%i\n" %
                                     (entry.mPredictionId, entry.mSbjctToken, entry.mSbjctStrand,
                                      old_f, old_t,
                                      f, t,
                                      first_res, last_res,
                                      mfirst_res, mlast_res,
                                      f, t))

                options.stderr.write(
                    "# %s\n" % str(map_sbjcts[entry.mSbjctToken]))
                options.stderr.write("# %s\n" % str(entry))
                options.stderr.flush()
                nerrors_boundaries += 1
                is_error = True

                # get extended boundaries for alignment later on
                while mfirst_res == 0 and first_res > 1:
                    first_res -= 1
                    mfirst_res = map_a2b.mapRowToCol(first_res)
                while mlast_res == 0 and last_res < map_a2b.getRowTo():
                    last_res += 1
                    mlast_res = map_a2b.mapRowToCol(last_res)

            # convert to genomic coordinates
            # convert negative strand coordinates
            if is_positive:
                new_f = mfirst_res - 1
                new_t = mlast_res
            else:
                new_f = mfirst_res
                new_t = mlast_res - 1

                new_f = map_a2b.getColTo() - new_f
                new_t = map_a2b.getColTo() - new_t

            # Now map the alignment.
            try:
                MapAlignment(entry, map_a2b)

            except ValueError:
                options.stderr.write("# alignment mapping not possible for prediction %i on %s %s:%i-%i -> %i-%i -> %i-%i -> %i-%i -> %i-%i -> %i-%i\n" %
                                     (entry.mPredictionId, entry.mSbjctToken, entry.mSbjctStrand,
                                      old_f, old_t,
                                      f, t,
                                      first_res, last_res,
                                      mfirst_res, mlast_res,
                                      new_f, new_t,
                                      entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))
                options.stderr.write(
                    "# %s\n" % str(map_sbjcts[entry.mSbjctToken]))
                options.stderr.flush()
                nerrors_map += 1
                is_error = True

            if new_f != entry.mSbjctGenomeFrom or new_t != entry.mSbjctGenomeTo:
                options.stderr.write("# mapping inconsistency for prediction %i on %s %s:%i-%i -> %i-%i -> %i-%i -> %i-%i -> %i-%i <> %i-%i\n" %
                                     (entry.mPredictionId, entry.mSbjctToken, entry.mSbjctStrand,
                                      old_f, old_t,
                                      f, t,
                                      first_res, last_res,
                                      mfirst_res, mlast_res,
                                      new_f, new_t,
                                      entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))

                nerrors_inconsistencies += 1
                is_error = True

        #######################################################################
        # Map via offsets
        if entry.mSbjctToken in breakpoints:

            old_f, old_t = entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo

            # convert to forward coordinates:
            if is_positive:
                f, t = old_f, old_t
            else:
                f, t = contig_sizes_old[
                    entry.mSbjctToken] - old_t, contig_sizes_old[entry.mSbjctToken] - old_f

            o1 = GetOffset(f,
                           breakpoints[entry.mSbjctToken],
                           endpoints[entry.mSbjctToken],
                           offsets[entry.mSbjctToken])
            o2 = GetOffset(t,
                           breakpoints[entry.mSbjctToken],
                           endpoints[entry.mSbjctToken],
                           offsets[entry.mSbjctToken])

            if o1 != o2:
                options.stderr.write("# break within gene %s\n" % str(entry))
                nerrors_map += 1
                is_error = True

            f += o1
            t += o2

            if not is_positive:
                f, t = contig_sizes_new[
                    entry.mSbjctToken] - t, contig_sizes_new[entry.mSbjctToken] - f

            entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo = f, t

            if entry.mSbjctGenomeFrom > entry.mSbjctGenomeTo:
                options.stderr.write(
                    "# mapping error: start after end %s\n" % str(entry))
                nerrors_map += 1
                is_error = True

        #######################################################################
        # do translation check, if genome is given
        if options.genome_file:
            genomic_sequence = Genomics.GetGenomicSequence(new_sbjct_token, entry.mSbjctStrand,
                                                           entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                                                           options.genome_file,
                                                           loglevel=0)

            map_peptide2translation, translation = Genomics.Alignment2PeptideAlignment(
                entry.mMapPeptide2Genome, entry.mQueryFrom, 0, genomic_sequence)

            if re.sub("X", "", translation) != re.sub("X", "", entry.mTranslation):
                options.stderr.write("# translation error for prediction %i on %s %s:%i-%i -> %i-%i <> %i-%i\n" %
                                     (entry.mPredictionId, entry.mSbjctToken, entry.mSbjctStrand,
                                      old_f, old_t,
                                      f, t,
                                      entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))
                if map_sbjcts:
                    options.stderr.write(
                        "# %s\n" % str(map_sbjcts[entry.mSbjctToken]))
                options.stderr.write(
                    "# old=%s\n# new=%s\n" % (entry.mTranslation, translation))
                options.stderr.write("# old=%s\n# new=%s\n" % (
                    entry.mAlignmentString, Genomics.Alignment2String(entry.mMapPeptide2Genome)))
                nerrors_translation += 1
                is_error = True

                if peptides and entry.mQueryToken in peptides:
                    naligned += 1

                    options.stdlog.write("# aligning: %s versus %s:%s: %i-%i\n" % (
                        entry.mQueryToken,
                        new_sbjct_token, entry.mSbjctStrand,
                        entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))

                    # do a quick reprediction
                    if entry.mQueryToken in peptides:
                        genomic_sequence = Genomics.GetGenomicSequence(new_sbjct_token, entry.mSbjctStrand,
                                                                       0, 0,
                                                                       genome_file=options.genome_pattern,
                                                                       loglevel=0)
                        predictor.mLogLevel = 0

                        result = predictor(entry.mQueryToken, peptides[entry.mQueryToken],
                                           entry.mSbjctToken, genomic_sequence,
                                           "--exhaustive --subopt FALSE --score '%s' " % str(
                                               80),
                                           new_f - 10, new_t + 10)
                        prediction_id = entry.mPredictionId
                        if result:
                            entry = result[0]
                            entry.mPredictionId = prediction_id
                            nerrors_realigned += 1
            else:
                if is_error:
                    nerrors_inconsequential += 1

        entry.mSbjctToken = new_sbjct_token

        # map query tokens
        if query_map_old2new:
            query_tokens = query_map_old2new[entry.mQueryToken]
        else:
            query_tokens = (entry.mQueryToken,)

        if options.skip_errors and is_error:
            continue

        for query_token in query_tokens:

            entry.mQueryToken = query_token

            prediction_id = entry.mPredictionId
            entry.mPredictionId = 0

            hid = Genomics.GetHID(str(entry))
            if hid in output:
                nduplicates += 1
                continue

            noutput += 1
            if options.renumber:
                prediction_id = noutput

            entry.mPredictionId = prediction_id

            options.stdout.write(str(entry) + "\n")
            options.stdout.flush()
            found_transcripts[entry.mQueryToken] = 1

    # write out found transcripts and genes
    nmissed_transcripts = 0
    missed_transcripts = []
    found_genes = {}
    if peptides:
        for x in peptides.keys():
            if x not in found_transcripts:
                nmissed_transcripts += 1
                missed_transcripts.append(x)
            else:
                found_genes[map_transcript2gene[x]] = 1

    missed_genes = {}
    nmissed_genes = 0
    if map_transcript2gene:

        for t in missed_transcripts:
            g = map_transcript2gene[t]
            if g not in found_genes:
                missed_genes[g] = 1
        nmissed_genes = len(missed_genes)

    if options.write_missed:
        outfile = open(options.write_missed, "w")
        for x in missed_transcripts:
            if x in unmapped:
                status = "unmapped"
            else:
                status = "mapped"
            outfile.write("%s\t%s\t%s\n" % ("transcript", x, status))
        for x in missed_genes:
            status = "unknown"
            outfile.write("%s\t%s\t%s\n" % ("gene", x, status))

        outfile.close()

    options.stdlog.write("# input=%i, output=%i, filtered=%i, nduplicates=%i, mapped=%i, errors=%i\n" % (
        ninput, noutput, nfiltered, nduplicates, nmapped, nerrors))
    options.stdlog.write("# errors: inconsequental=%i, boundaries=%i, mapping=%i, inconsistencies=%i, translation=%i, realigned=%i\n" % (
        nerrors_inconsequential, nerrors_boundaries, nerrors_map, nerrors_inconsistencies, nerrors_translation, nerrors_realigned))
    options.stdlog.write("# peptides: input=%i, found=%i, missed=%i, found_genes=%i, missed_genes=%i\n" % (
        len(peptides), len(found_transcripts), nmissed_transcripts, len(found_genes), nmissed_genes))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
