"""
=====================================
Rnaseq.py - Tools for RNAseq analysis
=====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Pipeline components - GO analysis

Tasks related to gene set GO analysis.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

"""

import CGAT.Experiment as E
import CGAT.CSV as CSV

import os
import shutil
import itertools
import glob
import collections

import numpy
import CGAT.GTF as GTF
import CGAT.BamTools as BamTools
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.rinterface as ri
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

import CGAT.Pipeline as P

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
P.getParameters(
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS

if os.path.exists("pipeline_conf.py"):
    E.info("reading additional configuration from pipeline_conf.py")
    execfile("pipeline_conf.py")

#############################################################
#############################################################
#############################################################
# UTR estimation
#############################################################
Utr = collections.namedtuple("Utr", "old new max status")


def buildUTRExtension(infile, outfile):
    '''build new utrs by building and fitting an HMM
    to reads upstream and downstream of known genes.

    Works on output of buildGeneLevelReadExtension.

    Known problems

    * the size of the extension is limited by the window size

    * introns within UTRs are ignored.

    * UTR extension might be underestimated for highly expressed genes
      as relative read counts drop off quickly, even though there is
      a good amount of reads still present in the UTR.

    The model

    The model is a three-state model::

        UTR --|--> notUTR --|--> otherTranscript --|
          ^---|      ^------|              ^-------|
                     ^-----------------------------|

    The chain starts in UTR and ends in notUTr or otherTranscript.

    The otherTranscript state models peaks of within the upstream/
    downstream region of a gene. These peaks might correspond to
    additional exons or unknown transcripts. Without this state,
    the UTR might be artificially extend to include these peaks.

    Emissions are modelled with beta distributions. These
    distributions permit both bimodal (UTR) and unimodal (notUTR)
    distribution of counts.

    Parameter estimation

    Parameters are derived from known UTRs within full length 
    territories.

    Transitions and emissions for the otherTranscript state
    are set heuristically:

       * low probabibily for remaining in state "otherTranscript".
           * these transcripts should be short.

       * emissions biased towards high counts - only strong signals
           will be considered.

       * these could be estimated from known UTRs, but I am worried
           UTR extensions then will be diluted.


    Alternatives

    The method could be improved.

        * base level resolution? 
            * longer chains result in more data and longer running times.
            * the averaging in windows smoothes the data, which might have
                a beneficial effect.

        * raw counts instead of scaled counts?
            * better model, as highly expressed genes should give more
                confident predictions.

    '''

    # the bin size , see gtf2table - can be cleaned from column names
    # or better set as options in .ini file
    binsize = 100
    territory_size = 15000

    # read gene coordinates
    geneinfos = {}
    for x in CSV.DictReader(IOTools.openFile(infile), dialect='excel-tab'):
        contig, strand, start, end = x['contig'], x[
            'strand'], int(x['start']), int(x['end'])
        geneinfos[x['gene_id']] = (contig, strand,
                                   start, end)

    infiles = [infile + ".readextension_upstream_sense.tsv.gz",
               infile + ".readextension_downstream_sense.tsv.gz"]

    outdir = os.path.join(PARAMS["exportdir"], "utr_extension")

    R('''suppressMessages(library(RColorBrewer))''')
    R('''suppressMessages(library(MASS))''')
    R('''suppressMessages(library(HiddenMarkov))''')

    # for upstream, downstream
    upstream_utrs, downstream_utrs = {}, {}

    all_genes = set()

    for filename, new_utrs in zip(infiles, (upstream_utrs, downstream_utrs)):

        E.info("processing %s" % filename)

        parts = os.path.basename(filename).split(".")

        data = R(
            '''data = read.table( gzfile( "%(filename)s"), header=TRUE, fill=TRUE, row.names=1)''' % locals() )

        ##########################################
        ##########################################
        ##########################################
        # estimation
        ##########################################
        # take only those with a 'complete' territory
        R('''d = data[-which( apply( data,1,function(x)any(is.na(x)))),]''')
        # save UTR
        R('''utrs = d$utr''' )
        # remove length and utr column
        R('''d = d[-c(1,2)]''')
        # remove those which are completely empty, logtransform or scale data
        # and export
        R('''lraw = log10( d[-which( apply(d,1,function(x)all(x==0))),] + 1 )''')

        utrs = R('''utrs = utrs[-which( apply(d,1,function(x)all(x==0)))]''' )
        scaled = R(
            '''lscaled = t(scale(t(lraw), center=FALSE, scale=apply(lraw,1,max) ))''' )
        exons = R('''lraw[,1]''')

        #######################################################
        #######################################################
        #######################################################
        # do the estimation:
        E.debug("estimation: utrs=%i, exons=%i, vals=%i, dim=%s" %
                (len(utrs), len(exons), len(scaled), R.dim(scaled)))
        # counts within and outside UTRs
        within_utr, outside_utr, otherTranscript = [], [], []
        # number of transitions between utrs
        transitions = numpy.zeros((3, 3), numpy.int)

        for x in xrange(len(utrs)):
            utr, exon = utrs[x], exons[x]

            # only consider genes with expression coverage
            # note: expression level is logscaled here, 10^1 = 10
            if exon < 0.1:
                continue

            # first row is column names, so x + 1
            values = list(scaled.rx(x + 1, True))

            utr_bins = utr // binsize
            nonutr_bins = (territory_size - utr) // binsize

            # build transition matrix
            transitions[0][0] += utr_bins
            transitions[0][1] += 1
            transitions[1][1] += nonutr_bins

            outside_utr.extend([x for x in values[utr_bins:] if x <= 0.5])

            # ignore exon and zero counts
            within_utr.extend([x for x in values[1:utr_bins] if x > 0.1])

            # add only high counts to otherTranscript emissions
            otherTranscript.extend([x for x in values[utr_bins:] if x > 0.5])

        # estimation for
        # 5% chance of transiting to otherTranscript
        transitions[1][2] = transitions[1][1] * 0.05
        # 10% chance of remaining in otherTranscript
        transitions[2][1] = 900
        transitions[2][2] = 100

        E.info("counting: (n,mean): within utr=%i,%f, outside utr=%i,%f, otherTranscript=%i,%f" %
               (len(within_utr), numpy.mean(within_utr),
                len(outside_utr), numpy.mean(outside_utr),
                len(otherTranscript), numpy.mean(otherTranscript)))

        ro.globalenv['transitions'] = R.matrix(transitions, nrow=3, ncol=3)
        R('''transitions = transitions / rowSums( transitions )''')
        ro.globalenv['within_utr'] = ro.FloatVector(within_utr[:10000])
        ro.globalenv['outside_utr'] = ro.FloatVector(outside_utr[:10000])
        ro.globalenv['otherTranscript'] = ro.FloatVector(
            otherTranscript[:10000])

        # estimate beta distribution parameters
        R('''doFit = function( data ) {
                   data[data == 0] = data[data == 0] + 0.001
                   data[data == 1] = data[data == 1] - 0.001
                   f = fitdistr( data, dbeta, list( shape1=0.5, shape2=0.5 ) )
                   return (f) }''' )

        fit_within_utr = R(
            '''fit_within_utr = suppressMessages(doFit( within_utr))''' )
        fit_outside_utr = R(
            '''fit_outside_utr = suppressMessages(doFit( outside_utr))''' )
        fit_other = R(
            '''fit_otherTranscript = suppressMessages(doFit( otherTranscript))''' )

        within_a, within_b = list(fit_within_utr.rx("estimate"))[0]
        outside_a, outside_b = list(fit_outside_utr.rx("estimate"))[0]
        other_a, other_b = list(fit_other.rx("estimate"))[0]

        E.info("beta estimates: within_utr=%f,%f outside=%f,%f, other=%f,%f" %
               (within_a, within_b, outside_a, outside_b, other_a, other_b))

        fn = ".".join((parts[0], parts[4], "fit", "png"))
        outfilename = os.path.join(outdir, fn)
        R.png(outfilename, height=1000, width=1000)

        R( '''par(mfrow=c(3,1))''' )
        R( '''x=seq(0,1,0.02)''')
        R( '''hist( within_utr, 50, col=rgb( 0,0,1,0.2) )''' )
        R( '''par(new=TRUE)''')
        R(
            '''plot( x, dbeta( x, fit_within_utr$estimate['shape1'], fit_within_utr$estimate['shape2']), type='l', col='blue')''')

        R( '''hist( outside_utr, 50, col=rgb( 1,0,0,0.2 ) )''' )
        R( '''par(new=TRUE)''')
        R(
            '''plot( x, dbeta( x, fit_outside_utr$estimate['shape1'], fit_outside_utr$estimate['shape2']), type='l', col='red')''')

        R( '''hist( otherTranscript, 50, col=rgb( 0,1,0,0.2 ) )''' )
        R( '''par(new=TRUE)''')
        R(
            '''plot( x, dbeta( x, fit_otherTranscript$estimate['shape1'], fit_otherTranscript$estimate['shape2']), type='l', col='green')''')
        R['dev.off']()

        #####################################################
        #####################################################
        #####################################################
        # build hmm
        # state 1 = UTR
        # state 2 = notUTR
        # state 3 = other transcript
        p = R('''betaparams = list( shape1=c(fit_within_utr$estimate['shape1'],
                                         fit_outside_utr$estimate['shape1'],
                                         fit_otherTranscript$estimate['shape1']),
                                shape2=c(fit_within_utr$estimate['shape2'],
                                         fit_outside_utr$estimate['shape2'],
                                         fit_otherTranscript$estimate['shape2'])) ''')
        R('''hmm = dthmm(NULL, transitions, c(1,0,0), "beta", betaparams )''' )

        E.info("fitting starts")
        #####################################################
        #####################################################
        #####################################################
        # fit to every sequence
        genes = R('''rownames(data)''')
        all_genes.update(set(genes))
        utrs = R('''data$utr''')
        exons = R('''data$exon''')
        nseqs = len(utrs)

        counter = E.Counter()

        for idx in xrange(len(utrs)):

            gene_id = genes[idx]

            old_utr = utrs[idx]

            if idx % 100 == 0:
                E.debug("processing gene %i/%i" % (idx, len(utrs)))

            counter.input += 1

            # do not predict if terminal exon not expressed
            if exons[idx] < 1:
                counter.skipped_notexpressed += 1
                new_utrs[gene_id] = Utr._make(
                    (old_utr, None, None, "notexpressed"))
                continue

            R('''obs = data[%i,][-c(1,2)]''' % (idx + 1) )
            # remove na
            obs = R('''obs = obs[!is.na(obs)]''' )
            if len(obs) <= 1 or max(obs) == 0:
                new_utrs[gene_id] = Utr._make(
                    (old_utr, None, None, "no observations"))
                continue

            # normalize
            R('''obs = obs / max(obs)''')
            # add small epsilon to 0 and 1 values
            R('''obs[obs==0] = obs[obs==0] + 0.001 ''')
            R('''obs[obs==1] = obs[obs==1] - 0.001 ''')
            R('''hmm$x = obs''')

            states = None
            try:
                states = list(R('''states = Viterbi( hmm )'''))
            except ri.RRuntimeError, msg:
                counter.skipped_error += 1
                new_utrs[gene_id] = Utr._make((old_utr, None, None, "fail"))
                continue

            max_utr = binsize * (len(states) - 1)

            # subtract 1 for last exon
            try:
                new_utr = binsize * (states.index(2) - 1)
                new_utrs[gene_id] = Utr._make(
                    (old_utr, new_utr, max_utr, "ok"))
                counter.success += 1
            except ValueError:
                new_utrs[gene_id] = Utr._make(
                    (old_utr, max_utr, max_utr, "max"))
                counter.maxutr += 1

    E.info("fitting: %s" % str(counter))

    outf = IOTools.openFile(outfile, "w")

    outf.write("\t".join(["gene_id", "contig", "strand", "status5", "status3"] +
                         ["%s_%s_%s" % (x, y, z) for x, y, z in itertools.product(
                             ("old", "new", "max"),
                             ("5utr", "3utr"),
                             ("length", "start", "end"))]) + "\n")

    def _write(coords, strand):

        start5, end5, start3, end3 = coords
        if strand == "-":
            start5, end5, start3, end3 = start3, end3, start5, end5

        if start5 is None:
            start5, end5, l5 = "", "", ""
        else:
            l5 = end5 - start5

        if start3 is None:
            start3, end3, l3 = "", "", ""
        else:
            l3 = end3 - start3

        return "\t".join(map(str, (l5, start5, end5,
                                   l3, start3, end3)))

    def _buildCoords(upstream, downstream, start, end):

        r = []
        if upstream:
            start5, end5 = start - upstream, start
        else:
            start5, end5 = None, None
        if downstream:
            start3, end3 = end, end + downstream
        else:
            start3, end3 = None, None

        return start5, end5, start3, end3

    for gene_id in all_genes:

        contig, strand, start, end = geneinfos[gene_id]

        outf.write("%s\t%s\t%s" % (gene_id, contig, strand))

        if gene_id in upstream_utrs:
            upstream = upstream_utrs[gene_id]
        else:
            upstream = Utr._make((None, None, None, "missing"))
        if gene_id in downstream_utrs:
            downstream = downstream_utrs[gene_id]
        else:
            downstream = Utr._make((None, None, None, "missing"))

        if strand == "-":
            upstream, downstream = downstream, upstream

        # output prediction status
        outf.write("\t%s\t%s" % (upstream.status, downstream.status))

        # build upstream/downstream coordinates
        old_coordinates = _buildCoords(
            upstream.old, downstream.old, start, end)
        new_coordinates = _buildCoords(
            upstream.new, downstream.new, start, end)

        # reconciled = take maximum extension of UTR
        max_coordinates = []
        # note that None counts as 0 in min/max.
        for i, d in enumerate(zip(old_coordinates, new_coordinates)):
            if i % 2 == 0:
                v = [z for z in d if z is not None]
                if v:
                    max_coordinates.append(min(v))
                else:
                    max_coordinates.append(None)
            else:
                max_coordinates.append(max(d))

        # convert to 5'/3' coordinates
        outf.write("\t%s\t%s\t%s\n" % (_write(old_coordinates, strand),
                                       _write(new_coordinates, strand),
                                       _write(max_coordinates, strand)))

    outf.close()

#########################################################################
#########################################################################
#########################################################################


def plotGeneLevelReadExtension(infile, outfile):
    '''plot reads extending beyond last exon.'''

    infiles = glob.glob(infile + ".*.tsv.gz")

    outdir = os.path.join(PARAMS["exportdir"], "utr_extension")

    R('''suppressMessages(library(RColorBrewer))''')
    R('''suppressMessages(library(MASS))''')
    R('''suppressMessages(library(HiddenMarkov))''')

    # the bin size , see gtf2table - could be cleaned from column names
    binsize = 100
    territory_size = 15000

    for filename in infiles:

        E.info("processing %s" % filename)

        parts = os.path.basename(filename).split(".")

        data = R(
            '''data = read.table( gzfile( "%(filename)s"), header=TRUE, fill=TRUE, row.names=1)''' % locals() )

        ##########################################
        ##########################################
        ##########################################
        # estimation
        ##########################################
        # take only those with a 'complete' territory
        R('''d = data[-which( apply( data,1,function(x)any(is.na(x)))),]''')
        # save UTR
        R('''utrs = d$utr''' )
        # remove length and utr column
        R('''d = d[-c(1,2)]''')
        # remove those which are completely empty, logtransform or scale data
        # and export
        R('''lraw = log10( d[-which( apply(d,1,function(x)all(x==0))),] + 1 )''')

        utrs = R('''utrs = utrs[-which( apply(d,1,function(x)all(x==0)))]''' )
        scaled = R(
            '''lscaled = t(scale(t(lraw), center=FALSE, scale=apply(lraw,1,max) ))''' )
        exons = R('''lraw[,1]''')

        if len(utrs) == 0:
            E.warn("no data for %s" % filename)
            continue

        #######################################################
        #######################################################
        #######################################################
        R('''myplot = function( reads, utrs, ... ) {
           oreads = t(data.matrix( reads )[order(utrs), ] )
           outrs = utrs[order(utrs)]
           image( 1:nrow(oreads), 1:ncol(oreads), oreads ,
                  xlab = "", ylab = "",
                  col=brewer.pal(9,"Greens"),
                  axes=FALSE)
           # axis(BELOW<-1, at=1:nrow(oreads), labels=rownames(oreads), cex.axis=0.7)
           par(new=TRUE)
           plot( outrs, 1:length(outrs), yaxs="i", xaxs="i", 
                 ylab="genes", xlab="len(utr) / bp", 
                 type="S", 
                 xlim=c(0,nrow(oreads)*%(binsize)i))
        }''' % locals())

        fn = ".".join((parts[0], parts[4], "raw", "png"))
        outfilename = os.path.join(outdir, fn)

        R.png(outfilename, height=2000, width=1000)
        R('''myplot( lraw, utrs )''' )
        R['dev.off']()

        # plot scaled data
        fn = ".".join((parts[0], parts[4], "scaled", "png"))
        outfilename = os.path.join(outdir, fn)

        R.png(outfilename, height=2000, width=1000)
        R('''myplot( lscaled, utrs )''' )
        R['dev.off']()

    P.touch(outfile)

#############################################################################
#############################################################################
#############################################################################


def filterAndMergeGTF(infile, outfile, remove_genes, merge=False):
    '''filter gtf file infile with gene ids in remove_genes
    and write to outfile.

    If *merge* is set, the resultant transcript models are merged by overlap.

    A summary file "<outfile>.summary" contains the number of transcripts that failed 
    various filters.

    A file "<outfile>.removed.tsv.gz" contains the filters that a transcript failed.
    '''

    counter = E.Counter()

    # write summary table
    outf = IOTools.openFile(outfile + ".removed.tsv.gz", "w")
    outf.write("gene_id\tnoverlap\tsection\n")
    for gene_id, r in remove_genes.iteritems():
        for s in r:
            counter[s] += 1
        outf.write("%s\t%i\t%s\n" % (gene_id,
                                     len(r),
                                     ",".join(r)))
    outf.close()

    # filter gtf file
    tmpfile = P.getTempFile(".")
    inf = GTF.iterator(IOTools.openFile(infile))

    genes_input, genes_output = set(), set()

    for gtf in inf:
        genes_input.add(gtf.gene_id)
        if gtf.gene_id in remove_genes:
            continue
        genes_output.add(gtf.gene_id)
        tmpfile.write("%s\n" % str(gtf))

    tmpfile.close()
    tmpfilename = tmpfile.name

    outf = IOTools.openFile(outfile + ".summary.tsv.gz", "w")
    outf.write("category\ttranscripts\n")
    for x, y in counter.iteritems():
        outf.write("%s\t%i\n" % (x, y))
    outf.write("input\t%i\n" % len(genes_input))
    outf.write("output\t%i\n" % len(genes_output))
    outf.write("removed\t%i\n" % (len(genes_input) - len(genes_output)))

    outf.close()

    # close-by exons need to be merged, otherwise
     # cuffdiff fails for those on "." strand

    if merge:
        statement = '''
        %(scriptsdir)s/gff_sort pos < %(tmpfilename)s
        | python %(scriptsdir)s/gtf2gtf.py
            --unset-genes="NONC%%06i"
            --log=%(outfile)s.log
        | python %(scriptsdir)s/gtf2gtf.py
            --merge-genes
            --log=%(outfile)s.log
        | python %(scriptsdir)s/gtf2gtf.py
            --merge-exons
            --merge-exons-distance=5
            --log=%(outfile)s.log
        | python %(scriptsdir)s/gtf2gtf.py
            --renumber-genes="NONC%%06i"
            --log=%(outfile)s.log
        | python %(scriptsdir)s/gtf2gtf.py
            --renumber-transcripts="NONC%%06i"
            --log=%(outfile)s.log
        | %(scriptsdir)s/gff_sort genepos 
        | gzip > %(outfile)s
        '''
    else:
        statement = '''
        %(scriptsdir)s/gff_sort pos < %(tmpfilename)s
        | gzip > %(outfile)s
        '''

    P.run()

    os.unlink(tmpfilename)


#############################################################
#############################################################
#############################################################
# running cufflinks
#############################################################
def runCufflinks(infiles, outfile):
    '''estimate expression levels in each set.
    '''

    gtffile, bamfile = infiles

    job_options = "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]

    track = os.path.basename(P.snip(gtffile, ".gtf.gz"))

    tmpfilename = P.getTempFilename(".")
    if os.path.exists(tmpfilename):
        os.unlink(tmpfilename)

    gtffile = os.path.abspath(gtffile)
    bamfile = os.path.abspath(bamfile)
    outfile = os.path.abspath(outfile)

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    # increase max-bundle-length to 4.5Mb due to Galnt-2 in mm9 with a 4.3Mb
    # intron.

    # AH: removed log messages about BAM record error
    # These cause logfiles to grow several Gigs and are
    # frequent for BAM files not created by tophat.
    # Error is:
    # BAM record error: found spliced alignment without XS attribute
    statement = '''mkdir %(tmpfilename)s;
    cd %(tmpfilename)s;
    cufflinks --label %(track)s
              --GTF <(gunzip < %(gtffile)s)
              --num-threads %(cufflinks_threads)i
              --frag-bias-correct %(bowtie_index_dir)s/%(genome)s.fa
              --library-type %(cufflinks_library_type)s
              %(cufflinks_options)s
              %(bamfile)s
    | grep -v 'BAM record error'
    >& %(outfile)s;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s.gtf.gz;
    gzip < isoforms.fpkm_tracking > %(outfile)s.fpkm_tracking.gz;
    gzip < genes.fpkm_tracking > %(outfile)s.genes_tracking.gz;
    '''

    P.run()

    shutil.rmtree(tmpfilename)

#########################################################################
#########################################################################
#########################################################################


def loadCufflinks(infile, outfile):
    '''load expression level measurements.'''

    track = P.snip(outfile, ".load")
    P.load(infile + ".genes_tracking.gz",
           outfile=track + "_genefpkm.load",
           options="--index=gene_id "
           "--ignore-column=tracking_id "
           "--ignore-column=class_code "
           "--ignore-column=nearest_ref_id")

    track = P.snip(outfile, ".load")
    P.load(infile + ".fpkm_tracking.gz",
           outfile=track + "_fpkm.load",
           options="--index=tracking_id "
           "--ignore-column=nearest_ref_id "
           "--rename-column=tracking_id:transcript_id")

    P.touch(outfile)


def runFeatureCounts(annotations_file,
                     bamfile,
                     outfile,
                     nthreads=4,
                     strand=2,
                     options=""):
    '''run feature counts on *annotations_file* with
    *bam_file*.
    
    If the bam-file is paired, paired-end counting
    is enabled and the bam file automatically sorted.
    '''

    # featureCounts cannot handle gzipped in or out files
    outfile = P.snip(outfile, ".gz")
    tmpdir = P.getTempDir()
    annotations_tmp = os.path.join(tmpdir,
                                   'geneset.gtf')
    bam_tmp = os.path.join(tmpdir,
                           bamfile)

    # -p -B specifies count fragments rather than reads, and both
    # reads must map to the feature
    # for legacy reasons look at feature_counts_paired
    if BamTools.isPaired(bamfile):
        # select paired end mode, additional options
        paired_options = "-p -B"
        # remove .bam extension
        bam_prefix = P.snip(bam_tmp, ".bam")
        # sort by read name
        paired_processing = \
            """samtools 
                sort -@ %(nthreads)i -n %(bamfile)s %(bam_prefix)s; 
            checkpoint; """ % locals()
        bamfile = bam_tmp 
    else:
        paired_options = ""
        paired_processing = ""

    job_options = "-pe dedicated %i" % nthreads

    # AH: what is the -b option doing?
    statement = '''mkdir %(tmpdir)s;
                   zcat %(annotations_file)s > %(annotations_tmp)s;
                   checkpoint;
                   %(paired_processing)s
                   featureCounts %(options)s
                                 -T %(nthreads)i
                                 -s %(strand)s
                                 -b
                                 -a %(annotations_tmp)s
                                 %(paired_options)s
                                 -o %(outfile)s
                                 %(bamfile)s
                    >& %(outfile)s.log;
                    checkpoint;
                    gzip -f %(outfile)s;
                    checkpoint;
                    rm -rf %(tmpdir)s
    '''

    P.run()
