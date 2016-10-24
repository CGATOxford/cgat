import sys
import re
import random
import collections
import itertools
import os
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Pipeline as P

'''
Filters a GTF to generate a geneset over which to count reads and look for
differential expression.

Options

-c --remove-contigs
Followed by a list of contigs to remove delimited by |
e.g. -c chrMT|chrALT removes all rows of the GTF from chrMT and chrALT

-k --keep-contigs
Followed by a list of contigs to keep delimited by |
e.g. -k chr1|chr12|chr3 only keeps rows of the GTF from chr1, chr12 and chr3

-f --filter
Any number of filters can be specified as -f or --filter followed by a string
laid out as follows:

To keep only lines with a specific value in a specified column:
-f COLUMN=VALUE
e.g.
-f 'source=protein_coding' will return only rows where source = protein_coding

To keep only lines without a specific value in a specified column:
-f COLUMN!=VALUE
e.g.
-f 'source!=protein_coding' returns only rows where source != protein_coding

To keep only lines with values from a list in a specified column:
-f COLUMN-in-VALUE1+VALUE2+VALUE3
e.g.
-f 'source-in-protein_coding+lincRNA' will return only rows where source =
   protein_coding or source = lincRNA

To keep only lines with values not in a list in a specified column
-f COLUMN-notin-VALUE1+VALUE2+VALUE3
e.g.
-f 'source-notin-protein_coding+lincRNA' will return only rows where source !=
   protein_coding and source != lincRNA

To keep only lines with values listed in a file
-f COLUMN-in_file-FILENAME
e.g.
-f 'source-in_file-a.tsv' will read the file a.tsv and return only rows of
the gtf where source is one of the values in this file

To keep only lines with values not listed in a file
-f COLUMN-notin_file-FILENAME
e.g.
-f 'source-notin_file-a.tsv' will read the file a.tsv and return only rows of
the gtf where source is not one of the values in this file

'''

# the old pipeline could also do the following
# but these were only used by cufflinks - do we still want them?
# merge exons seperated by small introns
# ignore very long introns

# remove repetitive RNA - we don't understand what this is?


def filterGTF(gtf, filterstring, tempout):
    if "!=" in filterstring:
        column, value = filterstring.split("!=")
        filtertype = "notequals"
    elif "=" in filterstring:
        column, value = filterstring.split("=")
        filtertype = "equals"
    elif "-in-" in filterstring:
        column, value = filterstring.split("-in-")
        value = value.split("+")
        filtertype = "in"
    elif "-notin-" in filterstring:
        column, value = filterstring.split("-notin-")
        value = value.split("+")
        filtertype = "notin"
    elif "-in_file-" in filterstring:
        column, value = filterstring.split("-in_file-")
        value = [line.strip() for line in IOTools.openFile(value)]
        filtertype = "in_file"
    elif "-notin_file-" in filterstring:
        column, value = filterstring.split("-notin_file-")
        value = [line.strip() for line in IOTools.openFile(value)]
        filtertype = "notin_file"

    gfile = IOTools.openFile(gtf)
    G = GTF.iterator(gfile)

    out = IOTools.openFile(tempout, "w")
    for item in G:
        D = item.asDict()
        D['contig'] = item.contig
        D['source'] = item.source
        D['feature'] = item.feature
        D['start'] = item.start
        D['end'] = item.end
        D['strand'] = item.strand
        D['frame'] = item.frame

        if filtertype == "equals":
            if D[column] == value:
                out.write("%s\n" % str(item))
        elif filtertype == "notequals":
            if D[column] != value:
                out.write("%s\n" % str(item))
        elif filtertype == "in" or 'in_file':
            if D[column] in value:
                out.write("%s\n" % str(item))
        elif filtertype == "notin" or 'notin_file':
            if D[column] not in value:
                out.write("%s\n" % str(item))
    out.close()
    gfile.close()


def removeNamedContigs(contigs):
    return """awk '$1 !~ /(%(contigs)s)$/' |""" % locals()


def keepOnlyNamedContigs(contigs):
    return """awk '$1 ~ /(%(contigs)s)$/' |""" % locals()


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-g", "--gtf", dest="gtf",
                      type="string",
                      help="path to input gtf")

    parser.add_option("-c", "--remove-contigs", dest="remove_contigs",
                      type="string",
                      help="contigs to remove, delimited by |")

    parser.add_option("-k", "--keep-contigs", dest="keep_contigs",
                      type="string",
                      help="""all contigs to keep, delimited by |.  Contigs
                      specified in --remove-contigs will still be removed""")

    parser.add_option("-o", "--outfile", dest="outfile",
                      type="string",
                      help="path to processed output gtf")

    parser.add_option("-f", "--filter", dest="filters",
                      type="string",
                      action="append",
                      help="""List of filters to apply to your GTF""")

    parser.set_defaults(
        remove_contigs=None,
        keep_contigs=None,
    )

    (options, args) = E.Start(parser)

    gtf = options.gtf

    if options.remove_contigs or options.keep_contigs:
        statement = 'zcat %s |' % gtf

        if options.remove_contigs:
            statement += removeNamedContigs(options.remove_contigs)

        if options.keep_contigs:
             statement += keepOnlyNamedContigs(options.keep_contigs)

        if options.outfile.endswith(".gz"):
            outfile = options.outfile
        else:
            outfile = options.outfile + ".gz"

        statement += "gzip > %s " % outfile

        os.system(statement)

    T1 = gtf
    for filterstring in options.filters:
        T2 = P.getTempFilename(".")
        T2 = T2 + ".gtf"
        filterGTF(T1, filterstring, T2)
        T1 = T2

    shutil.move(T2, options.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
