
import sys
import os
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF

'''
Filters a GTF to generate a geneset over which to count reads and look for
differential expression.

Example:
This would keep only lines from the gtf where gene_biotype is pseudogene or
protein_coding, contig is not chrM or chrUn and transcript_support is greater
than 3.

python makeGeneset.py -i refcoding.gtf.gz - o outfile.gtf.gz
-f "gene_biotype=pseudogene+protein_coding" -f "transcript_support-morethan-3"
-c "chrM|chrUn"

GTF "columns" are named, in order
0 contig
1 source
2 feature
3 start
4 end
5 ?
6 strand
7 frame
plus anything before an "=" sign in column 8, e.g. gene_source, exon_number

Options
-g --gtf
path to local gtf to filter

-p --gtfpath
path to online gtf to download and filter

-o --outfile
path to outfile

-c --remove-contigs
Followed by a list of contigs to remove delimited by |
e.g. -c chrMT|chrALT removes all rows of the GTF from chrMT and chrALT

-k --keep-contigs
Followed by a list of contigs to keep delimited by |
e.g. -k chr1|chr12|chr3 only keeps rows of the GTF from chr1, chr12 and chr3

-f --filter
Any number of filters can be specified as -f or --filter followed by a string
laid out as follows:

To keep only lines with values from a list in a specified column:
-f COLUMN=VALUE1+VALUE2+VALUE3
e.g.
-f 'source=protein_coding+lincRNA' will return only rows where source =
   protein_coding or source = lincRNA

To keep only lines with values not in a list in a specified column
-f COLUMN!=VALUE1+VALUE2+VALUE3
e.g.
-f 'source!=protein_coding+lincRNA' will return only rows where source !=
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

To keep only lines with a value greater than or less than a threshold in
a specified column
-f COLUMN-morethan-VALUE or COLUMN-lessthan-VALUE
e.g. -f "transcript_support-morethan-3" will return only rows with transcript
support greater than 3

'''

# the old pipeline could also do the following
# but these were only used by cufflinks - do we still want them?
# merge exons seperated by small introns
# ignore very long introns
# remove repetitive RNA - we don't understand what this is?


def getGTF(path):
    statement = "wget %(path)s" % locals()
    os.system(statement)


def filterGTF(gtf, filterstring, tempout):

    if "!=" in filterstring:
        column, value = filterstring.split("!=")
        value = value.split("+")
        filtertype = "notin"

    elif "=" in filterstring:
        column, value = filterstring.split("=")
        value = value.split("+")
        filtertype = "in"

    elif "-in_file-" in filterstring:
        column, value = filterstring.split("-in_file-")
        value = [line.strip() for line in IOTools.openFile(value)]
        filtertype = "in_file"

    elif "-notin_file-" in filterstring:
        column, value = filterstring.split("-notin_file-")
        value = [line.strip() for line in IOTools.openFile(value)]
        filtertype = "notin_file"

    elif "-morethan-" in filterstring:
        column, value = filterstring.split("-morethan-")
        value = float(value)
        filtertype = "morethan"

    elif "-lessthan-" in filterstring:
        column, value = filterstring.split("-lessthan-")
        value = float(value)
        filtertype = "lessthan"

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

        if filtertype == "in" or filtertype == 'in_file':
            if D[column] in value:
                out.write("%s\n" % str(item))
        elif filtertype == "notin" or filtertype == 'notin_file':
            if D[column] not in value:
                out.write("%s\n" % str(item))
        elif filtertype == "morethan":
            if float(D[column]) > value:
                out.write("%s\n" % str(item))
        elif filtertype == "lessthan":
            if float(D[column]) < value:
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

    parser.add_option("-p", "--gtfpath", dest="gtfpath",
                      type="string",
                      help="path to online gtf")

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

    if options.gtf:
        gtf = options.gtf
    elif options.gtfpath:
        getGTF(options.gtfpath)
        gtf = options.gtfpath.split("/")[-1]
    else:
        raise ValueError("Please provide a GTF or the path to an online GTF")

    if not options.outfile:
        raise ValueError("Please provide an output file name")

    d = 0
    if options.remove_contigs or options.keep_contigs:
        d += 1
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
    if options.filters:
        d += 1
        for filterstring in options.filters:
            T2 = P.getTempFilename(".")
            T2 = T2 + ".gtf"
            filterGTF(T1, filterstring, T2)
            T1 = T2

        shutil.move(T2, options.outfile)

    if d == 0:
        raise ValueError("No filters provided")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
