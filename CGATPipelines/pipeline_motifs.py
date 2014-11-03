"""==============
Motif pipeline
==============

:Author: Andreas Heger
:Release: $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

The motif pipeline runs a set of motif discovery and enrichment
analysis on a set of intervals.


   * Motif discovery using MEME
   * Motif detection using MAST
      * requires a set of known motifs

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The
pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>`
documentation).  The configuration file is organized by section and
the variables are documented within the file. In order to get a local
configuration file in the current directory, type::

    python <codedir>/pipeline_motifs.py config

The following sections and parameters probably should be changed from the default 
values:

.. todo::
   describe important parameters

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.


Input
-----

Intervals
+++++++++

Input are :term:`bed`-formatted files of intervals. Intervals should
be at least bed4 formatted, i.e., each interval should be labelled
(uniquely).

Reference motifs
++++++++++++++++

Reference motifs are described by fasta-formatted multiple alignments, see for
example Jaspar download. The motifs are build by running MEME on the file.

Reference motifs should end in the suffix ".motif.fasta", for example,
:file:`rxrvdr.motif.fasta`.

Requirements
------------

The pipeline requires the information from the following pipelines:

:doc:`pipeline_annotations`
   set the configuration variable :py:data:`annotations_database` and
   :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie_             |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+

Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz
   tar -xvzf pipeline_chipseq.tgz
   cd pipeline_chipseq
   python <srcdir>/pipeline_chipseq.py make full

.. note::

   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   bowtie
      bowtie_ - a read mapper

.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml

Code
====

"""
import sys
import shutil
import itertools
import re
import glob
import os
from ruffus import *
import logging as L
import sqlite3
import xml.etree.ElementTree

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

import PipelineMotifs as PipelineMotifs
import CGATPipelines.PipelineTracks as PipelineTracks

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'annotations_dir': ""})

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__")

###################################################################
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.Sample

TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(glob.glob("*.bed.gz"),
                                                         "(\S+).bed.gz")

TRACKS_BEDFILES = ["%s.bed.gz" % x for x in TRACKS]

###################################################################
###################################################################
###################################################################
# if conf.py exists: execute to change the above assignmentsn
if os.path.exists("pipeline_conf.py"):
    L.info("reading additional configuration from pipeline_conf.py")
    execfile("pipeline_conf.py")

###################################################################
###################################################################
###################################################################
#
###################################################################


def getAssociatedBAMFiles(track):
    '''return a list of BAM files associated with a track.

    By default, this method searches for ``track.bam`` 
    file in the current directory and returns an offset of 0.

    Associations can be defined in the .ini file in the section
    [bams]. For example, the following snippet associates track
    track1 with the bamfiles :file:`track1.bam` and :file:`track2.bam`::

       [bams]
       track1=track1.bam,track2.bam

    Glob expressions are permitted.

    Offsets are used to shift tags in ChIP experiments. Offsets
    need to be defined in the [offsets] sections. If no offsets
    are defined, the method returns a list of 0 offsets.

    Offsets need to be defined in the same order as the bam files::

       [offsets]
       track1=120,200

    returns a list of BAM files and offsets.

    Default tracks and offsets can be specified using a placeholder ``%``. The
    following will associate all tracks with the same bam file::

        [bams]
        %=all.bam


    '''
    fn = track.asFile()
    bamfiles = glob.glob("%s.bam" % fn)

    if bamfiles == []:
        if "bams_%s" % fn.lower() in PARAMS:
            for ff in P.asList(PARAMS["bams_%s" % fn.lower()]):
                bamfiles.extend(glob.glob(ff))
        else:
            for pattern, value in P.CONFIG.items("bams"):
                if "%" in pattern:
                    p = re.sub("%", "\S+", pattern)
                    if re.search(p, fn, re.IGNORECASE):
                        bamfiles.extend(glob.glob(value))

    offsets = []
    if "offsets_%s" % fn.lower() in PARAMS:
        offsets = map(int, P.asList(PARAMS["offsets_%s" % fn.lower()]))
    else:
        for pattern, value in P.CONFIG.items("offsets"):
            if "%" in pattern:
                p = re.sub("%", "\S+", pattern)
                if re.search(p, fn, re.IGNORECASE):
                    offsets.extend(map(int, value.split(",")))

    if offsets == []:
        offsets = [0] * len(bamfiles)

    if len(bamfiles) != len(offsets):
        raise ValueError("number of BAM files %s is not the same as number of offsets: %s" % (
            str(bamfiles), str(offsets)))

    return bamfiles, offsets

###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
# General preparation tasks
###################################################################

############################################################
############################################################
############################################################


@transform(TRACKS_BEDFILES,
           suffix(".bed.gz"),
           "_intervals.load")
def loadIntervals(infile, outfile):
    '''load intervals from :term:`bed` formatted files into database.
    '''

    bedfile = infile

    track = Sample(filename=P.snip(infile, ".bed.gz"))
    bamfiles, offsets = getAssociatedBAMFiles(track)
    control = ""

    if bamfiles:
        E.info("%s: associated bamfiles = %s" % (track, bamfiles))
    else:
        E.info("%s: no bamfiles associated" % (track))

    assert (len(bamfiles) == 1)
    bamfile = bamfiles[0]
    offset = offsets[0]

    tablename = P.toTable(outfile)

    statement = '''zcat %(bedfile)s
                | awk '{printf("%%s\\t%%i\\t%%i\\t%%i\\n", $1,$2,$3,++a)}'
                | python %(scriptsdir)s/bed2table.py 
                           --counter=peaks
                           --bam-file=%(bamfile)s
                           --offset=%(offset)i
                           --bed-header=contig,start,end,interval_id
                           %(control)s
                           --all-fields 
                           --log=%(outfile)s
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --add-index=contig,start
                       --add-index=interval_id
                       --table=%(tablename)s
                       --allow-empty-file 
                > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################


@follows(mkdir(os.path.join(PARAMS["exportdir"], "bed")))
@transform(TRACKS_BEDFILES,
           regex(r"(.*).bed.gz"),
           os.path.join(PARAMS["exportdir"], "bed", r"\1.bed.gz"))
def indexIntervals(infile, outfile):
    '''index intervals.
    '''
    statement = '''zcat %(infile)s | sort -k1,1 -k2,2n | bgzip > %(outfile)s; tabix -p bed %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################


@follows(mkdir(os.path.join(PARAMS["exportdir"], "peaks")))
@transform(loadIntervals,
           regex(r"(.*)_intervals.load"),
           os.path.join(PARAMS["exportdir"], "peaks", r"\1.peak.bed.gz"))
def exportPeakLocations(infile, outfile):
    '''export peak locations
    '''

    dbh = connect()
    outf = IOTools.openFile(outfile, "w")
    cc = dbh.cursor()
    table = P.toTable(infile)
    for x in cc.execute( """SELECT contig, peakcenter, peakcenter+1, interval_id, peakval 
                                   FROM %(table)s """ % locals() ):
        outf.write("\t".join(map(str, x)) + "\n")
    outf.close()


############################################################
############################################################
############################################################
@transform(loadIntervals,
           suffix("_intervals.load"),
           ".discovery.fasta")
def exportMotifDiscoverySequences(infile, outfile):
    '''export sequences for motif discovery.

    This method requires the _interval tables.

    For motif discovery, only the sequences with the highest S/N ratio are supplied.

    1. The top *motifs_proportion* intervals sorted by peakval
    2. Only a region +/- *motifs_halfwidth* around the peak 
    3. At least *motifs_min_sequences*. If there are not enough sequences
          to start with, all will be used.
    4. At most *motifs_max_size* sequences will be output.
    '''
    track = P.snip(infile, "_intervals.load")
    dbhandle = connect()

    p = P.substituteParameters(**locals())
    nseq = PipelineMotifs.writeSequencesForIntervals(track,
                                                     outfile,
                                                     dbhandle,
                                                     full=False,
                                                     masker=P.asList(
                                                         p['motifs_masker']),
                                                     halfwidth=int(
                                                         p["motifs_halfwidth"]),
                                                     maxsize=int(
                                                         p["motifs_max_size"]),
                                                     proportion=p[
                                                         "motifs_proportion"],
                                                     min_sequences=p[
                                                         "motifs_min_sequences"],
                                                     num_sequences=p[
                                                         "motifs_num_sequences"],
                                                     order=p['motifs_score'])

    if nseq == 0:
        E.warn("%s: no sequences - meme skipped" % outfile)
        P.touch(outfile)


############################################################
############################################################
############################################################
@follows(mkdir("motifs"))
@transform(TRACKS_BEDFILES,
           regex("(.*).bed.gz"),
           r"motifs/\1.foreground.fasta")
def exportMotifDetectionSequences(infile, outfile):
    '''export sequences for motif discovery.

    This method requires the _interval tables.
    '''
    PipelineMotifs.exportSequencesFromBedFile(infile, outfile,
                                              masker=PARAMS['motifs_masker'])


@follows(mkdir("motifs"))
@transform(TRACKS_BEDFILES,
           regex("(.*).bed.gz"),
           r"motifs/\1.control.fasta")
def exportMotifControlSequences(infile, outfile):
    '''for each interval, export the left and right 
    sequence segment of the same size.
    '''
    PipelineMotifs.exportSequencesFromBedFile(infile, outfile,
                                              masker=PARAMS['motifs_masker'],
                                              mode="leftright")


############################################################
############################################################
############################################################
@follows(mkdir("meme.dir"))
@transform(exportMotifDiscoverySequences,
           regex("(.*).discovery.fasta"),
           r"meme.dir/\1.meme")
def runMeme(infile, outfile):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''

    track = P.snip(infile, ".discovery.fasta")

    PipelineMotifs.runMEMEOnSequences(infile, outfile)

############################################################
############################################################
############################################################


@merge(runMeme, "meme_summary.load")
def loadMemeSummary(infiles, outfile):
    '''load information about motifs into database.'''

    outf = P.getTempFile(".")

    outf.write("track\n")

    for infile in infiles:
        if IOTools.isEmpty(infile):
            continue
        motif = P.snip(infile, ".meme")
        outf.write("%s\n" % motif)

    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)

###############################################################################
############################ Meme-Chip ########################################
###############################################################################


def suggestMotifDiscoveryForeground():
    '''output bed files for motif discovery.
    '''

    npeaks = [x.strip() for x in str(PARAMS["memechip_npeaks"]).split(",")]
    widths = [x.strip() for x in str(PARAMS["memechip_widths"]).split(",")]
    maskers = [x.strip() for x in str(PARAMS["memechip_maskers"]).split(",")]

    for infile in TRACKS_BEDFILES:
        track = P.snip(os.path.basename(infile), ".bed.gz")
        for n, w, masker in itertools.product(npeaks, widths, maskers):
            foreground = os.path.join(
                "discovery.dir", ".".join([track, n, w, masker, "foreground", "fasta"]))
            background = os.path.join(
                "discovery.dir", ".".join([track, n, w, masker, "background", "fasta"]))
            yield (track + "_intervals.load", foreground, int(n), int(w), masker)


def suggestMotifDiscoveryBackground():
    '''output bed files for motif discovery.
    '''

    npeaks = [x.strip() for x in str(PARAMS["memechip_npeaks"]).split(",")]
    widths = [x.strip() for x in str(PARAMS["memechip_widths"]).split(",")]
    maskers = [x.strip() for x in str(PARAMS["memechip_maskers"]).split(",")]

    for infile in TRACKS_BEDFILES:
        track = P.snip(os.path.basename(infile), ".bed.gz")
        for n, w, masker in itertools.product(npeaks, widths, maskers):
            background = os.path.join(
                "discovery.dir", ".".join([track, n, w, masker, "background", "fasta"]))
            yield (track + "_intervals.load", background, int(n), int(w), masker)


@follows(loadIntervals, mkdir("discovery.dir"))
@files(suggestMotifDiscoveryForeground)
def buildDiscoverySequences(infile, outfile, npeaks, width, masker):
    '''get the peak sequences, masking or not specificed in the ini file.
    '''

    track = P.snip(infile, "_intervals.load")
    dbhandle = connect()

    nseq = PipelineMotifs.writeSequencesForIntervals(track,
                                                     outfile,
                                                     dbhandle,
                                                     full=False,
                                                     masker=[masker],
                                                     halfwidth=width,
                                                     maxsize=int(
                                                         PARAMS["motifs_max_size"]),
                                                     proportion=None,
                                                     num_sequences=npeaks,
                                                     order='peakval')

    if nseq == 0:
        E.warn("%s: no sequences in foreground" % outfile)
        P.touch(outfile)


@follows(loadIntervals, mkdir("discovery.dir"))
@files(suggestMotifDiscoveryBackground)
def buildBackgroundSequences(infile, outfile, npeaks, width, masker):
    '''get the peak sequences, masking or not specificed in the ini file.
    '''

    track = P.snip(infile, "_intervals.load")
    dbhandle = connect()

    nseq = PipelineMotifs.writeSequencesForIntervals(track,
                                                     outfile,
                                                     dbhandle,
                                                     full=False,
                                                     masker=[masker],
                                                     halfwidth=width,
                                                     maxsize=int(
                                                         PARAMS["motifs_max_size"]),
                                                     proportion=None,
                                                     num_sequences=npeaks,
                                                     order='peakval',
                                                     shift="leftright")

    if nseq == 0:
        E.warn("%s: no sequences in background" % outfile_background)
        # P.touch( outfile )

############################################################
############################################################
############################################################


@transform(buildBackgroundSequences,
           suffix(".fasta"),
           ".markov")
def buildMemeBackgroundFiles(infile, outfile):
    '''prepare the meme background model'''
    statement = '''fasta-get-markov -m 2 %(infile)s  > %(outfile)s''' % locals()
    P.run()

############################################################
############################################################
############################################################


@follows(mkdir("memechip.dir"))
@merge(PARAMS["memechip_transfac_matrices"],
       "memechip.dir/transfac.filtered.dat")
def filterTransfac(infile, outfile):
    '''filter the transfac matrices, here for vertebrate'''

    statement = '''cat %(infile)s 
                | python %(scriptsdir)s/transfac2transfac.py 
                         --method=filter --filter-method=V
                         --log=%(outfile)s.log
                >  %(outfile)s
                '''
    P.run()

############################################################
############################################################
############################################################


@transform(filterTransfac,
           suffix(".dat"),
           ".meme")
def makeMemeMotifs(infile, outfile):
    '''convert transfac motifs to meme format'''

    statement = '''transfac2meme 
                    -use_acc 
                    -logodds
                    %(infile)s
                   > %(outfile)s
                '''
    P.run()


############################################################
############################################################
############################################################
@follows(mkdir("memechip.dir"))
@collate((buildDiscoverySequences, buildMemeBackgroundFiles),
         regex("discovery.dir/(.*).(foreground.fasta|background.markov)"),
         r"memechip.dir/\1.memechip")
def runMemeChip(infiles, outfile):

    background_markov, foreground_fasta = infiles

    assert foreground_fasta.endswith(".fasta")

    transfacMotifs = PARAMS["memechip_transfac_meme"]

    nmotifs = PARAMS["memechip_nmotifs"]

    ncpu = PARAMS["memechip_ncpu"]

    # job_options = "-pe mpi %(ncpu)i " % locals()
    # job_queue = "mpi.q"
    # '-meme-p %(ncpu)i'

    outdir = os.path.join(os.path.abspath(PARAMS["exportdir"]),
                          "memechip",
                          os.path.basename(outfile))

    # remove any existing output directory as otherwise meme will fail
    try:
        shutil.rmtree(outdir)
    except OSError:
        pass

    statement = '''meme-chip
                   -o %(outdir)s
                   -db %(transfacMotifs)s
                   -bfile %(background_markov)s
                   -ccut 0
                   -meme-mod zoops
                   -meme-minw %(memechip_minw)s
                   -meme-maxw %(memechip_maxw)s
                   -meme-nmotifs %(memechip_nmotifs)s
                   -meme-maxsize %(memechip_max_size)i
                   %(foreground_fasta)s
                   %(memechip_options)s
                   > %(outfile)s
                '''

    P.run()

############################################################
############################################################
############################################################


@merge(runMemeChip, "memechip_summary.load")
def loadMemeChipSummary(infiles, outfile):
    '''load information about motifs into database.'''

    outf = P.getTempFile(".")

    outf.write("track\tnpeaks\twidth\tmasking\tpath\n")

    for infile in infiles:
        if IOTools.isEmpty(infile):
            continue
        fn = P.snip(os.path.basename(infile), ".memechip")

        track, npeaks, width, masking = fn.split(".")
        outf.write(
            "\t".join(map(str, (track, npeaks, width, masking, fn))) + "\n")

    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)


############################################################
############################################################
############################################################
@transform(exportMotifDiscoverySequences,
           suffix(".fasta"),
           ".motifseq_stats.load")
def loadMotifSequenceComposition(infile, outfile):
    '''compute sequence composition of sequences used for ab-initio search.'''

    to_cluster = True

    tablename = P.toTable(outfile)

    statement = '''
    python %(scriptsdir)s/fasta2table.py 
        --section=na
        --log=%(outfile)s
    < %(infile)s
    | python %(scriptsdir)s/csv2db.py
        %(csv2db_options)s
        --table=%(tablename)s
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################


@merge("*.motif", "motif_info.load")
def loadMotifInformation(infiles, outfile):
    '''load information about motifs into database.'''

    outf = P.getTempFile(".")

    outf.write("motif\n")

    for infile in infiles:
        if IOTools.isEmpty(infile):
            continue
        motif = P.snip(infile, ".motif")
        outf.write("%s\n" % motif)

    outf.close()

    P.load(outf.name, outfile, "--allow-empty-file")

    os.unlink(outf.name)

############################################################
############################################################
############################################################
# run against database of known motifs
############################################################


@transform(runMeme, suffix(".meme"), ".tomtom")
def runTomTom(infile, outfile):
    '''compare ab-initio motifs against tomtom.'''
    PipelineMotifs.runTomTom(infile, outfile)


@transform(runTomTom, suffix(".tomtom"), "_tomtom.load")
def loadTomTom(infile, outfile):
    '''load tomtom results'''

    tablename = P.toTable(outfile)

    resultsdir = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), "tomtom", infile)
    xml_file = os.path.join(resultsdir, "tomtom.xml")

    if not os.path.exists(xml_file):
        E.warn("no tomtom output - skipped loading ")
        P.touch(outfile)
        return

    # get the motif name from the xml file

    tree = xml.etree.ElementTree.ElementTree()
    tree.parse(xml_file)
    motifs = tree.find("targets")
    name2alt = {}
    for motif in motifs.getiterator("motif"):
        name = motif.get("name")
        alt = motif.get("alt")
        name2alt[name] = alt

    tmpfile = P.getTempFile(".")

    # parse the text file
    for line in IOTools.openFile(infile):
        if line.startswith("#Query"):
            tmpfile.write(
                "target_name\tquery_id\ttarget_id\toptimal_offset\tpvalue\tevalue\tqvalue\tOverlap\tquery_consensus\ttarget_consensus\torientation\n")
            continue
        data = line[:-1].split("\t")
        target_name = name2alt[data[1]]
        tmpfile.write("%s\t%s" % (target_name, line))
    tmpfile.close()

    P.load(tmpfile.name, outfile)

    os.unlink(tmpfile.name)

############################################################
############################################################
############################################################


@files_re((exportMotifDetectionSequences, exportMotifControlSequences),
          "(\S+).control.fasta",
          [r"\1.control.fasta", r"\1.foreground.fasta",
           glob.glob("*.motif")],
          r"\1.mast.gz")
def runMast(infiles, outfile):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that
    all sequences are output and MAST curves can be computed. 

    10000 is a heuristic.
    '''
    PipelineMotifs.runMAST(infiles, outfile)

############################################################
############################################################
############################################################


@transform(runMast,
           suffix(".mast.gz"),
           "_mast.load")
def loadMast(infile, outfile):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''
    PipelineMotifs.loadMAST(infile, outfile)

############################################################
############################################################
############################################################


@follows(loadMotifInformation, mkdir(os.path.join(PARAMS["exportdir"], "motifs")))
@merge(loadMast, "motifs.export")
def exportMotifLocations(infiles, outfile):
    '''export motif locations. There will be a bed-file per motif.

    Overlapping motif matches in different tracks will be merged.
    '''

    dbh = connect()
    cc = dbh.cursor()

    motifs = [x[0]
              for x in cc.execute("SELECT motif FROM motif_info").fetchall()]

    for motif in motifs:

        tmpf = P.getTempFile(".")

        for infile in infiles:
            table = P.toTable(infile)
            track = P.snip(table, "_mast")
            for x in cc.execute( """SELECT contig, start, end, '%(track)s', evalue
                                   FROM %(table)s WHERE motif = '%(motif)s' AND start IS NOT NULL""" % locals() ):
                tmpf.write("\t".join(map(str, x)) + "\n")
        tmpf.close()

        outfile = os.path.join(
            PARAMS["exportdir"], "motifs", "%s.bed.gz" % motif)
        tmpfname = tmpf.name

        statement = '''mergeBed -i %(tmpfname)s -nms | gzip > %(outfile)s'''
        P.run()

        os.unlink(tmpf.name)


############################################################
############################################################
############################################################
# export section start
############################################################

############################################################
############################################################
############################################################
# export section end
############################################################


###################################################################
###################################################################
###################################################################
@follows(loadMemeChipSummary,
         loadIntervals)
def full():
    '''run the full pipeline.'''
    pass

###################################################################
###################################################################
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)

###################################################################
###################################################################
###################################################################


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)

###################################################################
###################################################################
###################################################################


@follows(mkdir("%s/bedfiles" % PARAMS["web_dir"]),
         update_report,
         )
def publish():
    '''publish files.'''
    # publish web pages

    P.publish_report()

    # publish additional data
    web_dir = PARAMS["web_dir"]
    project_id = P.getProjectId()

    # directory, files
    exportfiles = {
        "intervals": glob.glob(os.path.join(PARAMS["exportdir"], "bed", "*.bed.gz")) +
        glob.glob(os.path.join(PARAMS["exportdir"], "bed", "*.bed.gz.tbi")),
    }

    bams = []

    for targetdir, filenames in exportfiles.iteritems():
        if len(filenames) == 0:
            E.warn("no files for target '%s'" % targetdir)
        for src in filenames:
            dest = "%s/%s/%s" % (web_dir, targetdir, os.path.basename(src))
            if dest.endswith(".bam"):
                bams.append(dest)
            dest = os.path.abspath(dest)
            destdir = os.path.dirname(dest)
            if not os.path.exists(destdir):
                os.makedirs(destdir)

            if not os.path.exists(dest):
                E.debug("creating symlink from %s to %s" % (src, dest))
                os.symlink(os.path.abspath(src), dest)

    # output ucsc links
    for bam in bams:
        filename = os.path.basename(bam)
        track = P.snip(filename, ".bam")
        print """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/bamfiles/%(filename)s""" % locals()

if __name__ == "__main__":

    # print( "# tracks found: %s" % TRACKS_ALL )
    # print( "# tracks by experiment: %s" % TRACKS_EXPERIMENTS )
    # print( "# tracks by condition: %s" % TRACKS_CONDITIONS )
    # print( "# tracks by tissue: %s" % TRACKS_TISSUES )

    # fails for config command
    # check compatibility
    # assert PARAMS["genome"] == PARAMS_ANNOTATIONS["genome"]
    # sanity checks
    # assert len(TRACKS_ALL) > 0

    sys.exit(P.main(sys.argv))
