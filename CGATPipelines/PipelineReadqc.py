import os
import re
import glob
import cStringIO
import CGAT.Pipeline as P
import CGAT.IOTools as IOTools
import CGAT.CSV2DB as CSV2DB

###############################################################################
###############################################################################
###############################################################################
# Generate summaries from FastQC runs
###############################################################################


def FastqcSectionIterator(infile):
    data = []
    for line in infile:
        if line.startswith(">>END_MODULE"):
            yield name, status, header, data
        elif line.startswith(">>"):
            name, status = line[2:-1].split("\t")
            data = []
        elif line.startswith("#"):
            header = "\t".join([x for x in line[1:-1].split("\t") if x != ""])
        else:
            data.append(
                "\t".join([x for x in line[:-1].split("\t") if x != ""]))


def collectFastQCSections(infiles, section, datadir):
    '''iterate over all fastqc files and extract a particular section.'''
    results = []
    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        filename = os.path.join(datadir, track + "*_fastqc", "fastqc_data.txt")
        for fn in glob.glob(filename):
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                if name == section:
                    results.append((track, status, header, data))
    return results


def loadFastqc(filename):
    '''load FASTQC stats.'''

    for fn in glob.glob(filename):
        prefix = os.path.basename(os.path.dirname(fn))
        results = []

        for name, status, header, data in FastqcSectionIterator(
                IOTools.openFile(fn)):
            # do not collect basic stats, see loadFastQCSummary
            if name == "Basic Statistics":
                continue

            parser = CSV2DB.buildParser()
            (options, args) = parser.parse_args([])
            options.tablename = prefix + "_" + re.sub(" ", "_", name)
            options.allow_empty = True

            inf = cStringIO.StringIO("\n".join([header] + data) + "\n")
            CSV2DB.run(inf, options)
            results.append((name, status))

        # load status table
        parser = CSV2DB.buildParser()
        (options, args) = parser.parse_args([])
        options.tablename = prefix + "_status"
        options.allow_empty = True

        inf = cStringIO.StringIO(
            "\n".join(["name\tstatus"] +
                      ["\t".join(x) for x in results]) + "\n")
        CSV2DB.run(inf, options)


def buildFastQCSummaryStatus(infiles, outfile, datadir):
    '''load fastqc status summaries into a single table.'''

    outf = IOTools.openFile(outfile, "w")
    first = True
    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        filename = os.path.join(datadir,
                                track + "*_fastqc",
                                "fastqc_data.txt")

        for fn in glob.glob(filename):
            prefix = os.path.basename(os.path.dirname(fn))
            results = []

            names, stats = [], []
            for name, status, header, data in FastqcSectionIterator(
                    IOTools.openFile(fn)):
                stats.append(status)
                names.append(name)

            if first:
                outf.write("track\tfilename\t%s\n" % "\t".join(names))
                first = False

            outf.write("%s\t%s\t%s\n" %
                       (track, os.path.dirname(fn), "\t".join(stats)))
    outf.close()


def buildFastQCSummaryBasicStatistics(infiles, outfile, datadir):
    '''load fastqc summaries into a single table.'''

    data = collectFastQCSections(infiles, "Basic Statistics", datadir)

    outf = IOTools.openFile(outfile, "w")
    first = True
    for track, status, header, rows in data:
        rows = [x.split("\t") for x in rows]
        if first:
            headers = [row[0] for row in rows]
            outf.write("track\t%s\n" % "\t".join(headers))
            first = False
        outf.write("%s\t%s\n" % (track, "\t".join([row[1] for row in rows])))
    outf.close()
