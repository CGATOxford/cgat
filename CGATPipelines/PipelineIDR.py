import os
import re
import pysam
import random

import CGAT.Pipeline as P
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineTracks as PipelineTracks

###############################################################################
###############################################################################
###############################################################################
# Pre Process bamfiles
###############################################################################


def splitBam(infile, outfile_stub, params):
    pysam_in = pysam.Samfile(infile, "rb")
    n_outfiles = int(params[0])

    outfile_handles = []
    outfile_names = []

    # create list of upper bounds for intervals
    intervals = []
    lower = 0
    for i in range(n_outfiles):
        upper = lower + 1.0 / n_outfiles
        intervals.append(upper)
        # add an outfile handle to list of outfile handles
        outf = outfile_stub + "_" + str(i).zfill(2) + ".bam"
        outfile_handles.append(pysam.Samfile(outf, "wb", template=pysam_in))
        lower = upper

    # iterate through reads in samfile and write them to an outfile at random
    for read in pysam_in.fetch():
        r_num = random.random()
        for i in range(len(intervals)):
            if r_num < intervals[i]:
                outfile_handles[i].write(read)
                break
            else:
                continue

    # close outfiles
    for i in range(n_outfiles):
        outfile_handles[i].close()

    # index outfiles (pysam.index() throws an error)
    for split_sam in outfile_names:
        to_cluster = False
        statement = ("samtools index %(split_sam)s")
        P.run()


def filterBadLibraries(infiles, bad_samples):
    """
    Takes a list of infiles, removes those files that match any pattern in list
    of bad_samples, returns the filtered list of outfiles
    """
    bad_samples = [re.compile(x) for x in bad_samples]
    to_remove = []
    for inf in infiles:
        for regex in bad_samples:
            if regex.search(str(inf)):
                to_remove.append(inf)
    return to_remove


def mergeBams(infile_list, outfile):
    infile_list = " ".join(infile_list)
    out_stub = P.snip(outfile, ".bam")
    job_options = "-l mem_free=5G"
    statement = ("samtools merge - %(infile_list)s"
                 " | samtools sort - %(out_stub)s"
                 " 2>%(outfile)s.log;"
                 " checkpoint;"
                 " samtools index %(outfile)s"
                 " 2>%(outfile)s.bai.log")
    P.run()

##########################################################################
##########################################################################
##########################################################################
# Run Peak Calling For IDR
##########################################################################


def callIDRPeaks(infile,
                 outfile,
                 peak_caller,
                 control_option,
                 PARAMS_PEAKCALLER,
                 pseudoreplicate=False):

    # select peak caller to use
    if peak_caller == "macs2":
        caller = macs2IDRPeaks(control_option, PARAMS_PEAKCALLER)
    elif peak_caller == "spp":
        caller = sppIDRPeaks(control_option, PARAMS_PEAKCALLER)
    else:
        raise ValueError("Unrecognized peak-caller: %s" % peak_caller)

    caller.run(infile, outfile, pseudoreplicate)


class callerIDRPeaks(object):

    """
    Generic class for handling peak calling for IDR analysis
    """

    def __init__(self, control_option, PARAMS_PEAKCALLER):
        self.control_option = control_option
        self.PARAMS_PEAKCALLER = PARAMS_PEAKCALLER

    def getControlfile(self, track):
        """
        Return appropriate input file for a track.
        For pooled tracks (R0), will always return a pooled input file.
        If options set to pooled, will return a pooled input for all tracks.
        If options set to single, will return first input replicate (R1).
        If options set to matched, will return input with matching replicate.
        Otherwise will return ValueError
        """
        n = track.clone()
        n.data["attribute1"] = "input"  # is hardcoded into regex for ruffus tasks

        if n.replicate == "R0":
            # if track is pooled, then select pooled input
            control = os.path.basename("%s.bam" % n.asFile())
            control = os.path.join("./bamfiles_pooled", control)
        elif self.control_option == "pool":
            # if all controls are pooled, then select pooled input
            n.replicate = "R0"
            control = os.path.basename("%s.bam" % n.asFile())
            control = os.path.join("./bamfiles_pooled", control)
        elif self.control_option == "single":
            # if only one input file available, select it
            n.replicate = "R1"
            control = os.path.basename("%s.bam" % n.asFile())
            control = os.path.join("./bamfiles_filtered", control)
        elif self.control_option == "matching":
            # otherwise, select input with matching replicate
            control = os.path.basename("%s.bam" % n.asFile())
            control = os.path.join("bamfiles_filtered", control)
        else:
            raise ValueError("Unrecognised option_control %s"
                             " must be either pooled, single, or matching")

        if not os.path.exists(control):
            raise IOError("Control file is missing: %s" % control)

        return control

    def getRunStatement(self, infile, outfile,  controlfile):
        """
        Generate a specific run statement for each peakcaller class
        """
        return ""

    def postProcess(self, infile, outfile, controlfile):
        """
        Generate a specific post process statement for each peakcaller class
        that generates an outfile in the format required by IDR
        """
        return ""

    def run(self, infile, outfile, pseudoreplicate):
        """
        Gets gets appropriate input
        Generates a run statement
        Submits job
        Runs post processing steps for IDR
        """

        # get appropriate input
        Sample = PipelineTracks.AutoSample
        if pseudoreplicate:
            try:
                track = P.snip(infile, "_00.bam")
            except:
                track = P.snip(infile, "_01.bam")
        else:
            track = P.snip(infile, ".bam")
        controlfile = self.getControlfile(Sample(track))

        # following bugfix, check input file actually contains word 'input'
        assert re.search("input", os.path.basename(controlfile)), \
            "Input file doesn't contain 'input' in name: %s" % controlfile

        # run peakcalling
        statement = self.getRunStatement(infile, outfile, controlfile)
        # check run statement
        # print ("\nRun statement for sample %s :\n %s" % (infile, statement))
        job_options = "-l mem_free=10G"
        P.run()

        # post process peakcalling results
        ignore_pipe_errors = True
        statement = self.postProcess(infile, outfile, controlfile)
        # check post process statement
        # print ("\nPost-process statement for sample %s :"
        #        "\n %s" % (infile, statement))
        if statement:
            P.run()
        else:
            pass


class macs2IDRPeaks(callerIDRPeaks):

    """
    """

    def getRunStatement(self, infile, outfile, controlfile):
        """
        Generate a specific run statement for each peakcaller class
        """

        # generate outfile prefix
        dir_name = os.path.dirname(outfile)
        infile_stub = P.snip(os.path.basename(infile), ".bam")
        control_stub = P.snip(os.path.basename(controlfile), ".bam")
        outfile_stub = infile_stub + "_VS_" + control_stub
        outfile_stub = os.path.join(dir_name, outfile_stub)

        # build macs2 commandline statement
        statement = [("macs2 callpeak"
                      " --treatment %(infile)s"
                      " --control %(controlfile)s"
                      " --verbose=10")]

        # add additional parameters
        # currently the input read format has to be bam bc of ruffus regex
        statement.append("--format BAM")
        statement.append("--name %s" % outfile_stub)
        # require genome size, if it is not specified try to take from genome
        if not re.search("-g\s|--gsize",
                         self.PARAMS_PEAKCALLER["macs2_options_parameters"]):
            statement.append(
                "--gsize %s" % self.PARAMS_PEAKCALLER[
                    "macs2_options_genome_prefix"][:2])

        # set threshold for lax peak calling
        if self.PARAMS_PEAKCALLER["macs2_options_fdr"]:
            if self.PARAMS_PEAKCALLER["macs2_options_pvalue"]:
                raise Exception("Value specified for both macs2 options"
                                " -pvalue and -fdr please select one or"
                                " other option, but not both")
            else:
                threshold = "--qvalue " + \
                    str(self.PARAMS_PEAKCALLER["macs2_options_fdr"])
        elif self.PARAMS_PEAKCALLER["macs2_options_pvalue"]:
            threshold = "--pvalue=" + \
                str(self.PARAMS_PEAKCALLER["macs2_options_pvalue"])
        else:
            raise Exception("Must specify a value for either"
                            " macs2_options_pvalue or macs2_options_fdr,"
                            " but not both")
        statement.append(threshold)

        # deal with duplicate reads
        if self.PARAMS_PEAKCALLER["macs2_options_keep_duplicates"]:
            statement.append(
                "--keep-dup %s" % self.PARAMS_PEAKCALLER[
                    "macs2_options_keep_duplicates"])

        # add additional parameters
        statement.append(self.PARAMS_PEAKCALLER["macs2_options_parameters"])

        # write log information to sentinel file
        statement.append(">& %(outfile)s")

        statement = (" ".join(statement) % locals())

        return statement

    def postProcess(self, infile, outfile, controlfile):
        """
        Takes the narrowPeak files output by macs2.
        If macs2 given pvalue, then sorts by column 8 (-log10(pval))
        If macs2 given fdr, then sorts by column 9 (-log10(qval)).
        Generates a regionPeak File in the format required for IDR analysis
        N.B. IDR pipeline expects macs2 to output .encodePeak file, whereas it
        actually outputs a .narrowPeak file.
        """
        # generate outfile prefix
        dir_name = os.path.dirname(outfile)
        infile_stub = P.snip(os.path.basename(infile), ".bam")
        control_stub = P.snip(os.path.basename(controlfile), ".bam")
        outfile_stub = infile_stub + "_VS_" + control_stub
        outfile_stub = os.path.join(dir_name, outfile_stub)

        # set up sort statement
        if self.PARAMS_PEAKCALLER["macs2_options_pvalue"]:
            sort_column = "8nr,8nr"
        else:
            sort_column = "9nr,9nr"  # Check this!

        npeaks = self.PARAMS_PEAKCALLER["macs2_options_npeaks"]

        # format outfile as required for idr
        statement = [("sort"
                      "  -k %(sort_column)s"
                      "     %(outfile_stub)s_peaks.narrowPeak"
                      " | head -n %(npeaks)s"
                      " | gzip"
                      " > %(outfile_stub)s.regionPeak.gz")]
        # zip the original bedfile
        statement.append("; gzip %(outfile_stub)s_peaks.narrowPeak")
        # zip the excel file
        statement.append("; gzip %(outfile_stub)s_peaks.xls")
        # zip the summits file
        statement.append("; gzip %(outfile_stub)s_summits.bed")
        # create statement
        statement = ("".join(statement) % locals())

        return statement


class sppIDRPeaks(callerIDRPeaks):

    """
    Class for calling IDR peaks using spp.
    No postprocessing is run because spp outputs files in the format required
    for IDR
    """

    def getRunStatement(self, infile, outfile, controlfile):
        """
        Generate a specific run statement for each peakcaller class
        """
        # select location of the spp script to run
        if self.PARAMS_PEAKCALLER["spp_options_idr_script"] == "default":
            executable = P.which("run_spp.R")
        elif self.PARAMS_PEAKCALLER["spp_options_idr_script"] == "nodups":
            executable = P.which("run_spp_nodups.R")
        else:
            executable = self.PARAMS_PEAKCALLER["spp_options_idr_script"]
            try:
                os.path.exists(executable)
            except:
                raise IOError("SPP script not found: %s" % executable)

        # select the threshold for lax peak calling
        if self.PARAMS_PEAKCALLER["spp_options_npeaks"]:
            if self.PARAMS_PEAKCALLER["spp_options_fdr"]:
                raise Exception("Value specified for both SPP options"
                                " -npeaks and -fdr please select one or"
                                " other option, but not both")
            else:
                threshold = "-npeaks=" + \
                    str(self.PARAMS_PEAKCALLER["spp_options_npeaks"])
        elif self.PARAMS_PEAKCALLER["spp_options_fdr"]:
            threshold = "-fdr=" + \
                str(self.PARAMS_PEAKCALLER["spp_options_fdr"])
        else:
            raise Exception("Must specify a value for either"
                            " spp_options_npeaks or spp_options_fdr,"
                            " but not both")

        # build run statement for spp.
        # -savn is output.npeak.file (passed as NULL,
        #                             means filename based on infile)
        # -out is output.result.file
        # -odir defaults to os.path.dirname( infile )
        # -savn is save narrowpeak file
        # -savr is save regionpeak file
        #  (run_spp.R script throws an error if region peak is not output).
        statement = [("Rscript %(executable)s"
                      " -c=%(infile)s"
                      " -i=%(controlfile)s"
                      " %(threshold)s"
                      " -savn"
                      " -savr")]

        # add additional options
        statement.append(self.PARAMS_PEAKCALLER["spp_options_parameters"])

        # specify outfile
        statement.append(" -rf"
                         " -out=/stats/phantomPeakStatsReps.tab"
                         " >& %(outfile)s")

        statement = (" ".join(statement) % locals())

        return statement

##########################################################################
##########################################################################
##########################################################################
# Run IDR analysis
##########################################################################


def getIDRStatement(infile1,
                    infile2,
                    outfile,
                    overlap_ratio,
                    ranking_measure,
                    chr_table,
                    idr_wrapper):

    # get outfile stub
    inf1 = os.path.basename(infile1).split("_VS_")[0]
    inf2 = os.path.basename(infile2).split("_VS_")[0]
    out_prefix = os.path.join(os.path.dirname(outfile),
                              inf1 + "_vs_" + inf2)

    statement = ("python %(idr_wrapper)s"
                 "  --action=run"
                 "  --output-prefix=%(out_prefix)s"
                 "  --chromosome-table=%(chr_table)s"
                 "  --signal-value=%(ranking_measure)s"
                 "  --overlap-ratio=%(overlap_ratio)s"
                 "  %(infile1)s"
                 "  %(infile2)s"
                 "  >> %(outfile)s" % locals())

    return statement


def getIDRPlotStatement(infiles, outfile, idr_wrapper):
    """
    Receives list of infiles, the fist of which is a sentinel, the subsequent
    files are *uri.sav files output from run-batch-consistency.r script.
    Returns a run statement for batch-consistency-plot.r as it is wrapped in
    WrapperIDR.py
    """
    infile_prefixes = [P.snip(x, "-uri.sav") for x in infiles[1:]]
    infile_prefixes = " ".join(infile_prefixes)
    outfile_prefix = P.snip(outfile, ".pdf")

    statement = ("python %(idr_wrapper)s"
                 "  --action=plot"
                 "  --output-prefix=%(outfile_prefix)s"
                 "  %(infile_prefixes)s;"
                 " convert"
                 # "  -resize 125%%" why does this cause an error?
                 "  %(outfile_prefix)s.pdf"
                 "  %(outfile_prefix)s.png" % locals())

    return statement


##########################################################################
##########################################################################
##########################################################################
# Post Process PeakCalling
##########################################################################
def countPeaks(infiles, outf):
    """
    Count the number of peaks in each narrowPeak file
    """
    for infile in infiles:
        sample_id = os.path.basename(infile).split("_VS_")[0]
        tissue, condition, replicate = sample_id.split("-")
        experiment = tissue + "_" + condition
        n_peaks = str(len(IOTools.openFile(infile).readlines()))
        outf.write("\t".join([sample_id,
                              experiment,
                              tissue,
                              condition,
                              replicate,
                              n_peaks]) + "\n")
    outf.close()

##########################################################################
##########################################################################
##########################################################################
# Post Process IDR
##########################################################################


def findNPeaks(infiles, outfile, params):
    outf = IOTools.openFile(outfile, "w")
    outf.write("Tissue\t"
               "Condition\t"
               "Experiment\t"
               "idr_comp\t"
               "sample_1\t"
               "sample_2\t"
               "n_peaks\n")
    idr_threshold = float(params[0])

    # Hack: for only one infile, P.submit returns a string rather than a list
    if type(infiles) is str:
        infiles = [infiles, ]

    for inf in infiles:
        inf_name = P.snip(os.path.basename(inf), "-overlapped-peaks.txt")
        tissue = inf_name.split("-")[0]
        condition = inf_name.split("-")[1]
        experiment = "_".join([tissue, condition])
        sample1, sample2 = inf_name.split("_vs_")
        n_peaks = 0
        header = True
        for line in IOTools.openFile(inf):
            if header:
                header = False
                continue
            line = line.split()
            if float(line[10]) <= idr_threshold:
                n_peaks += 1
            else:
                continue
        outf.write(tissue + "\t"
                   + condition + "\t"
                   + experiment + "\t"
                   + inf_name + "\t"
                   + sample1 + "\t"
                   + sample2 + "\t"
                   + str(n_peaks) + "\n")

    outf.close()
