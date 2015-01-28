'''

Requirements:

* bedtools >= 2.21.0
* picardtools >= 1.106
* samtools >= 1.1
* MEDIPS >= 1.15.0

'''

import os
import re
import collections
import pandas
import math
import numpy
import numpy.ma as ma
import itertools
import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.BamTools as BamTools
import CGAT.IOTools as IOTools
import CGAT.Expression as Expression
import CGAT.Bed as Bed


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


def countTags(infile, outfile):
    '''count number of pairs in bed-file.'''

    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/bed2stats.py
    --per-contig
    --log=%(outfile)s.log
    >& %(outfile)s'''
    P.run()


def countReadsWithinWindows(bedfile,
                            windowfile,
                            outfile,
                            counting_method="midpoint",
                            memory_allocation="4G"):
    '''count reads given in *tagfile* within intervals in
    *windowfile*.

    Both files need to be :term:`bed` formatted.

    Counting is done using bedtools. The counting method
    can be 'midpoint' or 'nucleotide'.
    '''

    job_options = "-l mem_free=%s" % memory_allocation

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


def buildDMRStats(infiles, outfile, method, fdr_threshold=None):
    '''build dmr summary statistics.

    This method works from output files created by Expression.py
    (method="deseq" or method="edger") or runMEDIPS (method="medips")
    '''
    results = collections.defaultdict(lambda: collections.defaultdict(int))
    status = collections.defaultdict(lambda: collections.defaultdict(int))

    # deseq/edger
    f_significant = lambda x: x.significant == "1"
    f_up = lambda x: float(x.l2fold) > 0
    f_down = lambda x: float(x.l2fold) < 0
    f_fold2up = lambda x: float(x.l2fold) > 1
    f_fold2down = lambda x: float(x.l2fold) < -1
    f_key = lambda x: (x.treatment_name, x.control_name)
    f_status = lambda x: x.status

    outf = IOTools.openFile(outfile, "w")

    is_first = True
    for infile in infiles:

        xx = 0
        for line in IOTools.iterate(IOTools.openFile(infile)):
            key = f_key(line)

            r, s = results[key], status[key]
            r["tested"] += 1
            ss = f_status(line)
            s[ss] += 1

            if ss != "OK":
                continue

            is_significant = f_significant(line)
            up = f_up(line)
            down = f_down(line)
            fold2up = f_fold2up(line)
            fold2down = f_fold2down(line)
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

            if xx > 10000:
                break

        if is_first:
            is_first = False
            header1, header2 = set(), set()
            for r in results.values():
                header1.update(r.keys())
            for s in status.values():
                header2.update(s.keys())

            header = ["method", "treatment", "control"]
            header1 = list(sorted(header1))
            header2 = list(sorted(header2))

            outf.write("\t".join(header + header1 + header2) + "\n")

        for treatment, control in results.keys():
            key = (treatment, control)
            r = results[key]
            s = status[key]
            outf.write("%s\t%s\t%s\t" % (method, treatment, control))
            outf.write("\t".join([str(r[x]) for x in header1]) + "\t")
            outf.write("\t".join([str(s[x]) for x in header2]) + "\n")


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
    # define the two groups
    groups = sorted(set([x.group for x in design.values()]))

    # build a filtering statement
    groupA, groupB = groups

    def _buildMax(g, threshold):

        selected = [x for x, y in design.items() if y.group == g]
        if len(selected) > 1:
            return "max((%s)) < %f" % (
                ",".join(
                    ["int(r['%s'])" % x for x in selected]),
                threshold)
        elif len(selected) == 1:
            return "int(r['%s']) < %f" % (selected[0], threshold)
        else:
            raise ValueError("no groups found for 'g'" % g)

    def _buildSum(g, threshold):

        selected = [x for x, y in design.items() if y.group == g]
        if len(selected) > 1:
            return "sum((%s)) > %f" % (
                ",".join(
                    ["int(r['%s'])" % x for x in selected]),
                threshold)
        elif len(selected) == 1:
            return "int(r['%s']) > %f" % (selected[0], threshold)
        else:
            raise ValueError("no groups found for 'g'" % g)

    upper_levelA = _buildMax(groupA, max_per_sample)
    upper_levelB = _buildMax(groupB, max_per_sample)
    sum_levelA = _buildSum(groupA, sum_per_group)
    sum_levelB = _buildSum(groupB, sum_per_group)

    statement = '''
    zcat %(counts_file)s
    | python %(scriptsdir)s/csv_select.py
            --log=%(outfile)s.log
            "(%(upper_levelA)s and %(sum_levelB)s) or
             (%(upper_levelB)s and %(sum_levelA)s)"
    | python %(scriptsdir)s/runExpression.py
            --log=%(outfile)s.log
            --design-tsv-file=%(design_file)s
            --tags-tsv-file=-
            --method=mock
            --filter-min-counts-per-sample=0
    | gzip
    > %(outfile)s
    '''

    P.run()


def runDE(infiles, outfile, outdir,
          method="deseq",
          spike_file=None):
    '''run DESeq or EdgeR.

    The job is split into smaller sections. The order of the input
    data is randomized in order to avoid any biases due to chromosomes and
    break up local correlations.

    At the end, a new q-value is computed from all results.
    '''

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
    --output-filename-pattern=%(outdir)s/%%s
    --subdirs
    --output-regex-header="^test_id"
    "python %(scriptsdir)s/runExpression.py
              --method=%(method)s
              --tags-tsv-file=-
              --design-tsv-file=%(design_file)s
              --output-filename-pattern=%%DIR%%/%(prefix)s_
              --deseq-fit-type=%(deseq_fit_type)s
              --deseq-dispersion-method=%(deseq_dispersion_method)s
              --deseq-sharing-mode=%(deseq_sharing_mode)s
              --edger-dispersion=%(edger_dispersion)f
              --filter-min-counts-per-row=%(tags_filter_min_counts_per_row)i
              --filter-min-counts-per-sample=%(tags_filter_min_counts_per_sample)i
              --filter-percentile-rowsums=%(tags_filter_percentile_rowsums)i
              --log=%(outfile)s.log
              --fdr=%(edger_fdr)f"
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
                                  header=0,
                                  index_col=0)

    # normalize count column by total library size
    # have to explicitly convert data_frame to numpy
    # array with int64/float64 data type.  Otherwise
    # numpy.log will through an Attribute error (wrong
    # error to report) as it cannot handle python longs
    bed_frame = bed_frame.fillna(0.0)
    val_array = numpy.array(bed_frame.values, dtype=numpy.int64)
    geom_mean = geoMean(val_array)
    ratio_frame = bed_frame.apply(lambda x: x/geom_mean,
                                  axis=0)
    size_factors = ratio_frame.apply(numpy.median,
                                     axis=0)
    normalize_frame = bed_frame/size_factors
    # replace infs and -infs with Nas, then 0s
    normalize_frame.replace([numpy.inf, -numpy.inf], numpy.nan, inplace=True)
    normalize_frame = normalize_frame.fillna(0.0)
    normalize_frame.to_csv(outfile, sep="\t", index_label="interval")


def geoMean(array):
    '''
    Generate the geometric mean of a list or array,
    removing all zero-values but retaining total length
    '''
    if isinstance(array, pandas.core.frame.DataFrame):
        array = array.as_matrix()
    else:
        pass
    non_zero = ma.masked_values(array,
                                0)

    log_a = ma.log(non_zero)
    geom_mean = ma.exp(log_a.mean())

    return geom_mean


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

    foldchange = lambda x: math.log((x['3_y'] + 1.0)/(x['3_x'] + 1.0), 2)
    merge_frame[4] = merge_frame.apply(foldchange, axis=1)

    out_frame = merge_frame[[0, 1, 2, 4]]
    out_frame.to_csv(outfile,
                     sep="\t",
                     header=None,
                     index=None)


def runMEDIPSQC(infile, outfile):
    '''run MEDIPS QC targets.'''

    # note that the wrapper adds the filename
    # to the output filenames.
    job_options = "-l mem_free=10G"

    statement = """python %(scriptsdir)s/runMEDIPS.py
            --ucsc-genome=%(medips_genome)s
            --treatment=%(infile)s
            --toolset=saturation
            --toolset=coverage
            --toolset=enrichment
            --shift=%(medips_shift)s
            --extend=%(medips_extension)s
            --output-filename-pattern="medips.dir/%%s"
            --log=%(outfile)s.log
            | gzip
            > %(outfile)s
            """
    P.run()


def runMEDIPSDMR(design_file, outfile):
    '''run differential MEDIPS analysis according to
    designfile.
    '''
    job_options = "-l mem_free=30G"

    design = Expression.readDesignFile(design_file)

    # remove data tracks not needed
    design = [(x, y) for x, y in design.items() if y.include]

    # build groups
    groups = set([y.group for x, y in design])

    statements = []
    for pair1, pair2 in itertools.combinations(groups, 2):
        treatment = ["%s.bam" % x for x, y in design if y.group == pair1]
        control = ["%s.bam" % x for x, y in design if y.group == pair2]

        treatment = ",".join(treatment)
        control = ",".join(control)
        # outfile contains directory prefix
        statements.append(
            """python %(scriptsdir)s/runMEDIPS.py
            --ucsc-genome=%(medips_genome)s
            --treatment=%(treatment)s
            --control=%(control)s
            --toolset=dmr
            --shift=%(medips_shift)s
            --extend=%(medips_extension)s
            --window-size=%(medips_window_size)i
            --output-filename-pattern="%(outfile)s_%(pair1)s_vs_%(pair2)s_%%s"
            --fdr-threshold=%(medips_fdr)f
            --log=%(outfile)s.log
            > %(outfile)s.log2;
            checkpoint;
            zcat %(outfile)s_%(pair1)s_vs_%(pair2)s_data.tsv.gz
            | python %(scriptsdir)s/runMEDIPS.py
            --treatment=%(pair1)s
            --control=%(pair2)s
            --toolset=convert
            --fdr-threshold=%(medips_fdr)f
            --log=%(outfile)s.log
            | gzip
            > %(outfile)s
            """)

    P.run()


@P.cluster_runnable
def outputSpikeCounts(outfile, infile_name,
                      expression_nbins=None,
                      fold_nbins=None,
                      expression_bins=None,
                      fold_bins=None):

    df = pandas.read_csv(infile_name,
                         sep="\t",
                         index_col=0)

    E.debug("read %i rows and %i columns of data" % df.shape)

    if "edger" in outfile.lower():
        # edger: treatment_mean and control_mean do not exist
        # use supplied values directly.
        l10average = numpy.log(df['treatment_mean'])
        l2fold = numpy.log2(df['fold'])
    else:
        # use pseudocounts to compute fold changes
        treatment_mean = df['treatment_mean'] + 1
        control_mean = df['control_mean'] + 1
        # build log2 average values
        l10average = numpy.log((treatment_mean + control_mean) / 2)
        l2fold = numpy.log2(treatment_mean / control_mean)

    if expression_nbins is not None:
        mm = math.ceil(max(l10average))
        expression_bins = numpy.arange(0, mm, mm / expression_nbins)

    if fold_nbins is not None:
        mm = math.ceil(max(abs(min(l2fold)), abs(max(l2fold))))
        # ensure that range is centered on exact 0
        n = math.ceil(fold_nbins / 2.0)
        fold_bins = numpy.concatenate(
            (-numpy.arange(0, mm, mm / n)[:0:-1],
             numpy.arange(0, mm, mm / n)))

    # compute expression bins
    d2hist_counts, xedges, yedges = numpy.histogram2d(
        l10average, l2fold,
        bins=(expression_bins,
              fold_bins))

    dd = pandas.DataFrame(d2hist_counts)
    dd.index = list(xedges[:-1])
    dd.columns = list(yedges[:-1])
    dd.to_csv(IOTools.openFile(outfile, "w"),
              sep="\t")

    return df, d2hist_counts, xedges, yedges, l10average, l2fold


@P.cluster_runnable
def plotDETagStats(infiles, outfile):
    '''plot differential expression stats'''

    infile, composition_file = infiles
    Expression.plotDETagStats(
        infile, outfile,
        additional_file=composition_file,
        join_columns=("contig", "start", "end"),
        additional_columns=("CpG_density",
                            "length"))

    P.touch(outfile)


@P.cluster_runnable
def buildSpikeResults(infile, outfile):
    '''build matrices with results from spike-in and upload
    into database.

    The method will output several files:

    .spiked.gz: Number of intervals that have been spiked-in
               for each bin of expression and fold-change

    .power.gz: Global power analysis - aggregates over all
        ranges of fold-change and expression and outputs the
        power, the proportion of intervals overall that
        could be detected as differentially methylated.

        This is a table with the following columns:

        fdr - fdr threshold
        power - power level, number of intervals detectable
        intervals - number of intervals in observed data at given
                    level of fdr and power.
        intervals_percent - percentage of intervals in observed data
              at given level of fdr and power

    '''

    expression_nbins = 10
    fold_nbins = 10

    spikefile = P.snip(infile, '.tsv.gz') + '.spike.gz'

    if not os.path.exists(spikefile):
        E.warn('no spike data: %s' % spikefile)
        P.touch(outfile)
        return

    ########################################
    # output and load spiked results
    tmpfile_name = P.getTempFilename(shared=True)

    statement = '''zcat %(spikefile)s
    | grep -e "^spike" -e "^test_id"
    > %(tmpfile_name)s
    '''
    P.run()

    E.debug("outputting spiked counts")
    (spiked, spiked_d2hist_counts, xedges, yedges,
     spiked_l10average, spiked_l2fold) = \
        outputSpikeCounts(
            outfile=P.snip(outfile, ".power.gz") + ".spiked.gz",
            infile_name=tmpfile_name,
            expression_nbins=expression_nbins,
            fold_nbins=fold_nbins)

    ########################################
    # output and load unspiked results
    statement = '''zcat %(infile)s
    | grep -v -e "^spike"
    > %(tmpfile_name)s
    '''
    P.run()
    E.debug("outputting unspiked counts")

    (unspiked, unspiked_d2hist_counts, unspiked_xedges,
     unspiked_yedges, unspiked_l10average, unspiked_l2fold) = \
        outputSpikeCounts(
            outfile=P.snip(outfile, ".power.gz") + ".unspiked.gz",
            infile_name=tmpfile_name,
            expression_bins=xedges,
            fold_bins=yedges)

    E.debug("computing power")

    assert xedges.all() == unspiked_xedges.all()

    tmpfile = IOTools.openFile(tmpfile_name, "w")
    tmpfile.write("\t".join(
        ("expression",
         "fold",
         "fdr",
         "counts",
         "percent")) + "\n")

    fdr_thresholds = [0.01, 0.05] + list(numpy.arange(0.1, 1.0, 0.1))
    power_thresholds = numpy.arange(0.1, 1.1, 0.1)

    spiked_total = float(spiked_d2hist_counts.sum().sum())
    unspiked_total = float(unspiked_d2hist_counts.sum().sum())

    outf = IOTools.openFile(outfile, "w")
    outf.write("fdr\tpower\tintervals\tintervals_percent\n")

    # significant results
    for fdr in fdr_thresholds:
        take = spiked['qvalue'] < fdr

        # compute 2D histogram in spiked data below fdr threshold
        spiked_d2hist_fdr, xedges, yedges = \
            numpy.histogram2d(spiked_l10average[take],
                              spiked_l2fold[take],
                              bins=(xedges, yedges))

        # convert to percentage of spike-ins per bin
        spiked_d2hist_fdr_normed = spiked_d2hist_fdr / spiked_d2hist_counts
        spiked_d2hist_fdr_normed = numpy.nan_to_num(spiked_d2hist_fdr_normed)

        # set values without data to -1
        spiked_d2hist_fdr_normed[spiked_d2hist_counts == 0] = -1.0

        # output to table for database upload
        for x, y in itertools.product(range(len(xedges) - 1),
                                      range(len(yedges) - 1)):
            tmpfile.write("\t".join(map(
                str, (xedges[x], yedges[y],
                      fdr,
                      spiked_d2hist_fdr[x, y],
                      100.0 * spiked_d2hist_fdr_normed[x, y]))) + "\n")

        # take elements in spiked_hist_fdr above a certain threshold
        for power in power_thresholds:
            # select 2D bins at a given power level
            power_take = spiked_d2hist_fdr_normed >= power

            # select the counts in the unspiked data according
            # to this level
            power_counts = unspiked_d2hist_counts[power_take]

            outf.write("\t".join(map(
                str, (fdr, power,
                      power_counts.sum().sum(),
                      100.0 * power_counts.sum().sum()
                      / unspiked_total))) + "\n")

    tmpfile.close()
    outf.close()

    # upload into table
    method = P.snip(os.path.dirname(outfile), ".dir")
    tablename = P.toTable(
        P.snip(outfile, "power.gz") + method + ".spike.load")

    statement = '''cat %(tmpfile_name)s
    | python %(scriptsdir)s/csv2db.py
           --table=%(tablename)s
           --add-index=fdr
    > %(outfile)s.log'''

    P.run()
    os.unlink(tmpfile_name)
