import os
import re
import collections
import pandas
from math import log
import numpy as np
import random

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.BamTools as BamTools
import CGAT.IOTools as IOTools
import CGAT.Expression as Expression
import CGAT.Bed as Bed

#########################################################################
#########################################################################
#########################################################################


def convertReadsToIntervals(bamfile,
                            bedfile,
                            filtering_quality=None,
                            filtering_dedup=None,
                            filtering_dedup_method='picard'):
    '''convert reads in *bamfile* to *intervals*.

    This method converts read data into intervals for
    counting based methods.

    This method is not appropriated for RNA-Seq.

    Optional steps include:

    * deduplication - remove duplicate reads
    * quality score filtering - remove reads below a certain quality score.
    * paired ended data - merge pairs
    * paired ended data - filter by insert size

    '''
    track = P.snip(bedfile, ".bed.gz")

    is_paired = BamTools.isPaired(bamfile)
    current_file = bamfile
    tmpdir = P.getTempFilename()
    statement = ["mkdir %(tmpdir)s"]
    nfiles = 0

    if filtering_quality > 0:
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()
        statement.append('''samtools view
        -q %(filtering_quality)i -b
        %(current_file)s
        2>> %%(bedfile)s.log
        > %(next_file)s ''' % locals())

        nfiles += 1
        current_file = next_file

    if filtering_dedup is not None:
        # Picard's MarkDuplicates requries an explicit bam file.
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()

        if filtering_dedup_method == 'samtools':
            statement.append('''samtools rmdup - - ''')

        elif filtering_dedup_method == 'picard':
            statement.append('''MarkDuplicates
            INPUT=%(current_file)s
            OUTPUT=%(next_file)s
            ASSUME_SORTED=TRUE
            METRICS_FILE=%(bedfile)s.duplicate_metrics
            REMOVE_DUPLICATES=TRUE
            VALIDATION_STRINGENCY=SILENT
            2>> %%(bedfile)s.log ''' % locals())

        nfiles += 1
        current_file = next_file

    if is_paired:
        statement.append('''cat %(current_file)s
            | python %(scriptsdir)s/bam2bed.py
              --merge-pairs
              --min-insert-size=%(filtering_min_insert_size)i
              --max-insert-size=%(filtering_max_insert_size)i
              --log=%(bedfile)s.log
              -
            | python %(scriptsdir)s/bed2bed.py
              --method=sanitize-genome
              --genome-file=%(genome_dir)s/%(genome)s
              --log=%(bedfile)s.log
            | cut -f 1,2,3,4
            | sort -k1,1 -k2,2n
            | bgzip > %(bedfile)s''')
    else:
        statement.append('''cat %(current_file)s
            | python %(scriptsdir)s/bam2bed.py
              --log=%(bedfile)s.log
              -
            | python %(scriptsdir)s/bed2bed.py
              --method=sanitize-genome
              --genome-file=%(genome_dir)s/%(genome)s
              --log=%(bedfile)s.log
            | cut -f 1,2,3,4
            | sort -k1,1 -k2,2n
            | bgzip > %(bedfile)s''')

    statement.append("tabix -p bed %(bedfile)s")
    statement.append("rm -rf %(tmpdir)s")
    statement = " ; ".join(statement)
    P.run()

    os.unlink(tmpdir)

#########################################################################
#########################################################################
#########################################################################


def countReadsWithinWindows(bedfile,
                            windowfile,
                            outfile,
                            counting_method="midpoint"):
    '''count reads given in *tagfile* within intervals in
    *windowfile*.

    Both files need to be :term:`bed` formatted.

    Counting is done using bedtools. The counting method
    can be 'midpoint' or 'nucleotide'.
    '''
    job_options = "-l mem_free=4G"

    if counting_method == "midpoint":
        f = '''| awk '{a = $2+($3-$2)/2;
        printf("%s\\t%i\\t%i\\n", $1, a, a+1)}' '''
    elif counting_method == "nucleotide":
        f = ""
    else:
        raise ValueError("unknown counting method: %s" % counting_method)

    statement = '''
    zcat %(bedfile)s
    %(f)s
    | coverageBed -a stdin -b %(windowfile)s -split
    | sort -k1,1 -k2,2n
    | gzip
    > %(outfile)s
    '''

    P.run()


#########################################################################
#########################################################################
#########################################################################
def aggregateWindowsReadCounts(infiles,
                               outfile,
                               regex="(.*)\..*"):
    '''aggregate several results from coverageBed
    into a single file.

    *regex* is used to extract the track name from the filename.
    The default removes any suffix.

    coverageBed outputs the following columns:
    1 Contig
    2 Start
    3 Stop
    4 Name
    5 The number of features in A that overlapped (by at least one
      base pair) the B interval.
    6 The number of bases in B that had non-zero coverage from features in A.
    7 The length of the entry in B.
    8 The fraction of bases in B that had non-zero coverage from
      features in A.

    For bed: use column 5
    For bed6: use column 7
    For bed12: use column 13

    Windows without any counts will not be output.
    '''

    # get bed format
    bed_columns = Bed.getNumColumns(infiles[0])
    # +1 as awk is 1-based
    column = bed_columns - 4 + 1

    src = " ".join(['''<( zcat %s |
              awk '{printf("%%s:%%i-%%i\\t%%i\\n", $1,$2,$3,$%s );}' ) ''' %
                    (x, column) for x in infiles])
    tmpfile = P.getTempFilename(".")
    statement = '''paste %(src)s > %(tmpfile)s'''
    P.run()

    # build track names
    tracks = [re.search(regex, os.path.basename(x)).groups()[0]
              for x in infiles]

    outf = IOTools.openFile(outfile, "w")
    outf.write("interval_id\t%s\n" % "\t".join(tracks))

    for line in open(tmpfile, "r"):
        data = line[:-1].split("\t")
        genes = list(set([data[x] for x in range(0, len(data), 2)]))
        values = [int(data[x]) for x in range(1, len(data), 2)]
        if sum(values) == 0:
            continue
        assert len(genes) == 1, \
            "paste command failed, wrong number of genes per line: '%s'" % line
        outf.write("%s\t%s\n" % (genes[0], "\t".join(map(str, values))))

    outf.close()

    os.unlink(tmpfile)


#########################################################################
#########################################################################
#########################################################################
def buildDMRStats(infile, outfile, method):
    '''build dmr summary statistics.
    '''
    results = collections.defaultdict(lambda: collections.defaultdict(int))

    status = collections.defaultdict(lambda: collections.defaultdict(int))
    x = 0
    for line in IOTools.iterate(IOTools.openFile(infile)):
        key = (line.treatment_name, line.control_name)
        r, s = results[key], status[key]
        r["tested"] += 1
        s[line.status] += 1

        is_significant = line.significant == "1"
        up = float(line.l2fold) > 0
        down = float(line.l2fold) < 0
        fold2up = float(line.l2fold) > 1
        fold2down = float(line.l2fold) < -1
        fold2 = fold2up or fold2down

        if up:
            r["up"] += 1
        if down:
            r["down"] += 1
        if fold2up:
            r["l2fold_up"] += 1
        if fold2down:
            r["l2fold_down"] += 1

        if is_significant:
            r["significant"] += 1
            if up:
                r["significant_up"] += 1
            if down:
                r["significant_down"] += 1
            if fold2:
                r["fold2"] += 1
            if fold2up:
                r["significant_l2fold_up"] += 1
            if fold2down:
                r["significant_l2fold_down"] += 1

    header1, header2 = set(), set()
    for r in results.values():
        header1.update(r.keys())
    for s in status.values():
        header2.update(s.keys())

    header = ["method", "treatment", "control"]
    header1 = list(sorted(header1))
    header2 = list(sorted(header2))

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join(header + header1 + header2) + "\n")

    for treatment, control in results.keys():
        key = (treatment, control)
        r = results[key]
        s = status[key]
        outf.write("%s\t%s\t%s\t" % (method, treatment, control))
        outf.write("\t".join([str(r[x]) for x in header1]) + "\t")
        outf.write("\t".join([str(s[x]) for x in header2]) + "\n")

#########################################################################
#########################################################################
#########################################################################


def buildFDRStats(infile, outfile, method):
    '''compute number of windows called at different FDR.
    '''

    data = pandas.read_csv(IOTools.openFile(infile), sep="\t", index_col=0)

    assert data['treatment_name'][0] == data['treatment_name'][-1]
    assert data['control_name'][0] == data['control_name'][-1]

    treatment_name, control_name = data[
        'treatment_name'][0], data['control_name'][0]

    key = (treatment_name, control_name)
    fdrs = (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

    for fdr in fdrs:
        print "fdr"
        take = data['qvalue'] <= fdr

        significant = sum(take)
        print significant


def outputAllWindows(infile, outfile):
    '''output all Windows as a bed file with the l2fold change
    as a score.
    '''
    outf = IOTools.openFile(outfile, "w")
    for line in IOTools.iterate(IOTools.openFile(infile)):
        outf.write("\t".join(
            (line.contig, line.start, line.end,
             "%6.4f" % float(line.l2fold))) + "\n")

    outf.close()


def outputRegionsOfInterest(infiles, outfile,
                            max_per_sample=10, sum_per_group=40):
    '''output windows according to various filters.

    The output is a mock analysis similar to a differential expression
    result.

    '''
    job_options = "-l mem_free=64G"

    design_file, counts_file = infiles

    design = Expression.readDesignFile(design_file)

    # remove tracks not included in the design
    design = dict([(x, y) for x, y in design.items() if y.include])
    des_i = design_items()
    # define the two groups
    groups = sorted(set([x.group for x in design.values()]))

    # build a filtering statement
    groupA, groupB = groups
    upper_levelA = "max( (%s) ) < %f" % (
        ",".join(
            ["int(r['%s'])" % x for x, y in des_i if y.group == groupA]),
        max_per_sample)

    sum_levelA = "sum( (%s) ) > %f" % (
        ",".join(
            ["int(r['%s'])" % x for x, y in des_i if y.group == groupB]),
        sum_per_group)

    upper_levelB = "max( (%s) ) < %f" % (
        ",".join(
            ["int(r['%s'])" % x for x, y in des_i if y.group == groupB]),
        max_per_sample)

    sum_levelB = "sum( (%s) ) > %f" % (
        ",".join(
            ["int(r['%s'])" % x for x, y in des_i if y.group == groupA]),
        sum_per_group)

    statement = '''
    zcat %(counts_file)s
    | python %(scriptsdir)s/csv_select.py
            --log=%(outfile)s.log
            "(%(upper_levelA)s and %(sum_levelA)s) or
                            (%(upper_levelB)s and %(sum_levelB)s)"
    | python %(scriptsdir)s/runExpression.py
            --log=%(outfile)s.log
            --filename-design=%(design_file)s
            --filename-tags=-
            --method=mock
            --filter-min-counts-per-sample=0
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################


def runDE(infiles, outfile, outdir,
          method="deseq",
          spike_file=None):
    '''run DESeq or EdgeR.

    The job is split into smaller sections. The order of the input
    data is randomized in order to avoid any biases due to chromosomes and
    break up local correlations.

    At the end, a new q-value is computed from all results.
    '''

    to_cluster = True

    design_file, counts_file = infiles

    if spike_file is None:
        statement = "zcat %(counts_file)s"
    else:
        statement = '''python %(scriptsdir)s/combine_tables.py
                           --missing-value=0
                           --cat=filename
                           --log=%(outfile)s.log
                           %(counts_file)s %(spike_file)s
              | python %(scriptsdir)s/csv_cut.py
                           --remove filename
                           --log=%(outfile)s.log
        '''

    prefix = os.path.basename(outfile)

    # the post-processing strips away the warning,
    # renames the qvalue column to old_qvalue
    # and adds a new qvalue column after recomputing
    # over all windows.
    statement += '''
              | perl %(scriptsdir)s/randomize_lines.pl -h
              | %(cmd-farm)s
                  --input-header
                  --output-header
                  --split-at-lines=200000
                  --cluster-options="-l mem_free=8G"
                  --log=%(outfile)s.log
                  --output-pattern=%(outdir)s/%%s
                  --subdirs
              "python %(scriptsdir)s/runExpression.py
              --method=%(method)s
              --filename-tags=-
              --filename-design=%(design_file)s
              --output-filename-pattern=%%DIR%%/%(prefix)s_
              --deseq-fit-type=%(deseq_fit_type)s
              --deseq-dispersion-method=%(deseq_dispersion_method)s
              --deseq-sharing-mode=%(deseq_sharing_mode)s
              --filter-min-counts-per-row=%(tags_filter_min_counts_per_row)i
              --filter-min-counts-per-sample=%(tags_filter_min_counts_per_sample)i
              --filter-percentile-rowsums=%(tags_filter_percentile_rowsums)i
              --log=%(outfile)s.log
              --fdr=%(edger_fdr)f"
              | grep -v "warnings"
              | perl %(scriptsdir)s/regtail.pl ^test_id
              | perl -p -e "s/qvalue/old_qvalue/"
              | python %(scriptsdir)s/table2table.py
              --log=%(outfile)s.log
              --method=fdr
              --column=pvalue
              --fdr-method=BH
              --fdr-add-column=qvalue
              | gzip
              > %(outfile)s '''

    P.run()


def normalizeBed(infile, outfile):
    '''
    Normalize counts in a bed file to total library size.
    Return as a bedGraph format.  Written as a function
    to use P.submit to send to cluster.
    '''

    bed_frame = pandas.read_table(infile,
                                  sep="\t",
                                  compression="gzip",
                                  header=None,
                                  index_col=0)

    # normalize count column by total library size

    med = np.median(bed_frame[4])
    normalize = lambda x: (x)/(float(med) + 1.0) + random.randint(0, 1)
    bed_frame[7] = bed_frame[4].apply(normalize)

    bed_frame.to_csv(outfile, sep="\t",
                     header=None,
                     columns=[1, 2, 7])


def enrichmentVsInput(infile, outfile):
    '''
    Calculate the fold enrichment of the test data
    vs. the input data
    '''

    test_frame = pandas.read_table(infile[1],
                                   sep="\t",
                                   compression="gzip",
                                   header=None,
                                   index_col=None)

    input_frame = pandas.read_table(infile[0],
                                    sep="\t",
                                    compression="gzip",
                                    header=None,
                                    index_col=None)
    merge_frame = pandas.merge(test_frame,
                               input_frame,
                               how='left',
                               left_on=[0, 1, 2],
                               right_on=[0, 1, 2])

    foldchange = lambda x: log((x['3_y'] + 1.0)/(x['3_x'] + 1.0), 2)
    merge_frame[4] = merge_frame.apply(foldchange, axis=1)

    out_frame = merge_frame[[0, 1, 2, 4]]
    out_frame.to_csv(outfile,
                     sep="\t",
                     header=None,
                     index=None)
