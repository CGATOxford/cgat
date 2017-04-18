'''add_transcripts_gtf.py - manipulate transcript models
============================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Genomics Genesets GTF Manipulation

Purpose
-------

This script take a gene set in a :term:`gtf` format, adds
additional transcript models, and outputs a new gene set in
:term:`gtf` format to stdout. Note that this does not work from stdin
as the input has to be parsed twice.

The type of additional transcript_mode is chosen by the ``--method``
command line option.

methods:

  - skip-exons
  - incomplete
  - 3prime

Type::

    cgat add_transcripts_gtf --help

for command line options.

Command line Options
--------------------

'''
import copy
import sys
import random
import collections

import numpy as np

import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):

    if not argv:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])
    parser.add_option(
        "--method", dest="method",
        type="choice",
        choices=("skip_exons", "incomplete", "3prime"),
        help="Method by which new transcript models are introduced "
        "[%default].")

    parser.add_option(
        "--fraction", dest="fraction",
        type="float",
        help="Fraction new transcript models to add, relative to the"
        "number of transcript models in the input GTF. This can be"
        "greater than 1 "
        "[%default].")

    parser.add_option(
        "--min-3prime", dest="--min_3prime",
        type="int",
        help="minimum addition 5 prime bases [%default].")

    parser.add_option(
        "--max-3prime", dest="--max_3prime",
        type="int",
        help="minimum addition 5 prime bases [%default].")

    parser.add_option(
        "--infile-geneset-gtf", dest="infile_gtf",
        type="string",
        help="path for infile gtf "
        "[%default].")

    parser.add_option(
        "--outfile-geneset-tsv", dest="outfile_tsv",
        type="string",
        help="path for tsv with transcript info "
        "[%default].")


    parser.set_defaults(
        method='skip_exons',
        fraction=0,
        infile_gtf=None,
        outfile_tsv=None,
        min_3prime=100,
        max_3prime=1000
    )

    (options, args) = E.Start(parser, argv=argv)

    ninput, noutput, nshort = 0, 0, 0

    if options.method is None:
        raise ValueError("please specify a --method")

    if options.fraction <= 0:
        raise ValueError("fraction must be >0")

    if options.infile_gtf is None:
        raise ValueError("please specify --infile-geneset-gtf")

    if options.outfile_tsv is None:
        raise ValueError("please specify --outfile-geneset-tsv")

    info_out = IOTools.openFile(options.outfile_tsv, "w")

    # need to add check that additional 3' end fits into contig length
    if options.method == '3prime':

        info_out.write("%s\t%s\n" % ("transcript_id", "added_bases"))

        gtfs = GTF.iterator(IOTools.openFile(options.infile_gtf))

        transcripts = set()

        for gtf_lines in GTF.transcript_iterator(gtfs, strict=0):

            if gtf_lines[0].feature == 'exon':
                ninput += 1
                transcript_id = gtf_lines[0].transcript_id
                transcripts.add(transcript_id)

        skip_transcripts = list(
                np.random.choice(list(transcripts),
                                 int(ninput * options.fraction),
                                 replace=False))

        gtfs = GTF.iterator(IOTools.openFile(options.infile_gtf))
        for gtf_lines in GTF.transcript_iterator(gtfs, strict=0):
                
            if gtf_lines[0].feature == 'exon':

                transcript_id = gtf_lines[0].transcript_id

                # write out every transcript_id already in gtf
                noutput += 1

                for gtf in gtf_lines:
                    options.stdout.write("%s\n" % str(gtf))

                # if transcript has been selected, create a new model
                # for each time it has been selected
                transcript_id = gtf_lines[0].transcript_id

                info_out.write("%s\t%i\n" % (transcript_id, 0))

                if transcript_id not in skip_transcripts:
                    continue

                bases_added = int(np.random.uniform(
                    options.min_3prime, options.max_3prime))

                new_transcript_id = "%s_longer3prime_%i" % (transcript_id, bases_added)

                info_out.write("%s\t%i\n" % (new_transcript_id, bases_added))

                noutput += 1

                for gtf in gtf_lines[:-1]:
                    options.stdout.write("%s\n" % str(gtf))

                gtf2 = copy.copy(gtf_lines[-1])
                gtf2.setAttribute('transcript_id', new_transcript_id)
                gtf2.end += bases_added

                options.stdout.write("%s\n" % str(gtf2))

    elif options.method == 'skip_exons' or options.method == 'incomplete':

        # set minimum require number of exons and write tsv header
        if options.method == 'skip_exons':
            min_exons = 2
            info_out.write("transcript_id\texons\tskipped_exon\n")

        elif options.method == 'incomplete':
            min_exons = 4
            info_out.write("transcript_id\texons\tstart\tend\n")

        transcript2exons = collections.defaultdict(set)

        gtfs = GTF.iterator(IOTools.openFile(options.infile_gtf))

        for gtf_lines in GTF.transcript_iterator(gtfs, strict=0):

            if gtf_lines[0].feature == 'exon':
                ninput += 1

                transcript_id = gtf_lines[0].transcript_id

                if len(gtf_lines)>min_exons:
                    transcript2exons[transcript_id] = range(1, len(gtf_lines)+1)
                else:
                    nshort += 1

        # Make list with each transcript repeated for each available exon
        if options.method == 'skip_exons':

            flat_transcripts = [t_id for t_id in transcript2exons.keys()
                                for exon in transcript2exons[t_id][1:-1]]

        elif options.method == 'incomplete':
            flat_transcripts = [t_id for t_id in transcript2exons.keys()
                                for exon in transcript2exons[t_id][:-2]]

        try:
            skip_transcripts = list(
                np.random.choice(flat_transcripts, int(ninput * options.fraction),
                                     replace=False))
        except ValueError as e:
            raise ValueError(
                'cannot ask for more %s additional '
                'models than there are exons available: available exons=%i, '
                'additional_models=%i.   numpy error: %s' % (
                    options.method, len(flat_transcripts),
                    int(ninput * options.fraction), e))

        gtfs = GTF.iterator(IOTools.openFile(options.infile_gtf))

        for gtf_lines in GTF.transcript_iterator(gtfs, strict=0):

            if gtf_lines[0].feature == 'exon':
                
                transcript_id = gtf_lines[0].transcript_id
                exons = len(gtf_lines)

                if options.method == 'skip_exons':
                    info_out.write("%s\t%i\t0\n" % (transcript_id, exons))
                elif options.method == 'incomplete':
                    info_out.write("%s\t%i\t1\t%i\n" % (transcript_id, exons, exons))
            
                # write out every transcript_id already in gtf
                noutput += 1

                for gtf in gtf_lines:
                    options.stdout.write("%s\n" % str(gtf))

                # if transcript has been selected, create a new model
                # for each time it has been selected
                transcript_id = gtf_lines[0].transcript_id

                if transcript_id not in skip_transcripts:
                    continue

                if options.method == 'skip_exons':
                    remove_exons = np.random.choice(
                        transcript2exons[transcript_id][1:-1], 
                        skip_transcripts.count(transcript_id),
                        replace=False)

                    for skipped_exon in remove_exons:

                        new_transcript_id = "%s_skipped_%i_%i" % (
                            transcript_id, skipped_exon,
                            max(transcript2exons[transcript_id]))

                        info_out.write("%s\t%i\t%i\n" % (
                            new_transcript_id, exons, skipped_exon))

                        noutput += 1


                        for ix, gtf in enumerate(gtf_lines):

                            # we want to skip one exon
                            if (ix + 1) != skipped_exon:

                                gtf2 = copy.copy(gtf)
                                gtf2.setAttribute('transcript_id', new_transcript_id)

                                options.stdout.write("%s\n" % str(gtf2))
                
                

                elif options.method == 'incomplete':

                    start_exons = np.random.choice(
                        transcript2exons[transcript_id][:-2],
                        skip_transcripts.count(transcript_id),
                        replace=False)

                    for start in start_exons:

                        end = np.random.choice(
                            list(transcript2exons[transcript_id])[start+1:])

                        new_transcript_id = "%s_incomplete_%i_%i_%i" % (
                            transcript_id, start, end, max(transcript2exons[transcript_id]))

                        info_out.write("%s\t%i\t%i\t%i\n" % (
                            new_transcript_id, exons, start, end))

                        noutput += 1

                        for ix, gtf in enumerate(gtf_lines):

                            # we only want to include the exons between start and end
                            if (ix + 1) >= start:
                                if (ix+1) <= end:

                                    gtf2 = copy.copy(gtf)
                                    gtf2.setAttribute('transcript_id',
                                                      new_transcript_id)

                                    options.stdout.write("%s\n" % str(gtf2))

                                # if ix>end skip remaining exons
                                else:
                                    break

    E.info("ninput=%i, noutput=%i, nshort=%i" %
           (ninput, noutput, nshort))
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
