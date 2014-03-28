'''classes and utility functions for pipeline_genomeassembly.py

Different assembly tools will use different inputs. Some can take
fasta files whereas others will take fastq and in either case can be
paired-end (in the same or different files) or single end

'''

import os
import collections
import shutil

import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGAT.FastaIterator as FastaIterator
import CGAT.Fastq as Fastq
import glob
import CGAT.Metaphlan as Metaphlan
import numpy as np


class Format:

    '''
    class for assessing formats
    '''

    def fileFormat(self, infile):
        '''
        return the file format for the short read data
        can be one of
        fasta
        fastq
        fasta.gz
        fastq.gz
        fasta.1.gz
        '''
        possible_formats = [
            "fa", "fasta", "fastq", "fasta.gz", "fastq.gz", "fasta.1.gz", "fastq.1.gz"]
        format = None
        for f in possible_formats:
            if infile.endswith(f):
                format = f
        assert format, "file %s is not of correct format" % infile
        return format


class PairedData(Format):

    '''
    class for assessing paired end data
    '''

    def __init__(self):

        self.format = None
        self.paired = False
        self.paired_interleaved = False
        self.paired_separate = False

    def getTrack(self, infile):
        '''
        return the track for the file
        '''
        return P.snip(os.path.basename(infile), ".%s" % self.fileFormat(infile))

    def getFormat(self, infile):
        self.format = self.fileFormat(infile)
        return self.format

    def checkPairedFile(self, infile):
        '''
        return if the paired-end data file
        attached to the input file for first in pair reads
        '''
        format = self.getFormat(infile)
        track = self.getTrack(infile)
        read2 = track + ".fastq.2.gz"
        assert len(glob.glob(read2)) > 0, "cannot find %s file for read 2 in the pair for %s" % (
            read2, infile)
        return ("separate", read2)

    def checkPairs(self, infile, ntries=10):
        '''
        Function to check if pairs exist interleaved
        within a file. If not it will check for separate
        files for the read pair
        '''
        format = self.getFormat(infile)
        paired = False
        inf = IOTools.openFile(infile)
        pairs = set()
        if format in ["fasta", "fasta.gz"]:
            iterator = FastaIterator.iterator
        elif format in ["fastq", "fastq.gz"]:
            iterator = Fastq.iterate
        elif format in ["fasta.1.gz", "fastq.1.gz"]:
            return self.checkPairedFile(infile)

        c = 0
        for record in iterator(inf):
            if record.quals:
                # make sure there are not other "/" in the sequence name
                seq_id = record.identifier.split("/")
            else:
                seq_id = record.title.split("/")

            assert len(
                seq_id) < 3, "cannot deal with this sequence name %s" % seq_id
            seq_id = seq_id[0]
            if seq_id not in pairs:
                pairs.add(seq_id)
                if c >= ntries:
                    break
                paired = False
                return paired
            else:
                print "found pair for %s" % (seq_id)
                paired = "interleaved"
                break
        return paired

#############################
# pooling reads across
# conditions
#############################
def pool_reads(infiles, outfile):
    '''
    pool raw reads across conditions
    '''
    # check first file for format and pairedness
    paired = PairedData().checkPairs(infiles[0])
    format = PairedData().getFormat(infiles[0])
    infs = " ".join(infiles)

    if not paired:
        statement = '''zcat %(infs)s | gzip > %(outfile)s''' % locals()
    else:
        infs2 = " ".join(
            [P.snip(x, format) + format.replace("1", "2") for x in infiles])
        outf2 = outfile.replace(format, format.replace("1", "2"))
        statement = '''zcat %(infs)s | gzip > %(outfile)s;
                       zcat %(infs2)s | gzip > %(outf2)s''' % locals()
    return statement

#############################
# filtering contigs by length
#############################
def filterContigs(infile, outfile, length):
    '''
    filter contigs by length
    '''
    outf = open(outfile, "w")
    for fasta in FastaIterator.iterate(IOTools.openFile(infile)):
        seq_length = len(fasta.sequence)
        if seq_length < length:
            continue
        outf.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
    outf.close()

############################
# function for performing
# claculation of stats
############################
def contig_to_stats(contigs_file, stats_file, params):
    '''
    calculate descriptive stats for a set
    of contigs / scaffolds
    '''

    PARAMS = params

    if PARAMS["filter"]:
        f = PARAMS["filter"]
    else:
        f = 0

    # iterate over the contigs/scaffolds and return stats
    number_of_scaffolds = 0

    N = PARAMS["scaffold_n"]
    scaffold_lengths = []

    inf = open(contigs_file)
    for record in FastaIterator.iterate(inf):
        scaffold_length = len(list(record.sequence))
        if scaffold_length >= f:
            number_of_scaffolds += 1
            scaffold_lengths.append(scaffold_length)

    # mean, median and max contig/scaffold lengths
    mean_length = np.mean(scaffold_lengths)
    median_length = np.median(scaffold_lengths)
    max_length = max(scaffold_lengths)

    # iterate over contigs/scaffolds sorted by longest
    # and caculate the NX
    index = 0
    cum_length = 0
    total_length = sum(scaffold_lengths)
    for length in sorted(scaffold_lengths, reverse=True):
        while cum_length <= total_length * (float(N) / 100):
            index += 1
            cum_length += length

    # output the results
    outf = open(stats_file, "w")
    outf.write(
        "nscaffolds\tscaffold_length\tN%i\tmedian_length\tmean_length\tmax_length\n" % N)
    outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (number_of_scaffolds, total_length, sorted(
        scaffold_lengths, reverse=True)[index], str(median_length), str(mean_length), str(max_length)))

###############################
###############################
###############################
def build_scaffold_lengths(contigs_file, outfile, params):
    '''
    output the distribution of scaffold lengths
    '''
    inf = open(contigs_file)
    outf = open(outfile, "w")
    outf.write("scaffold_name\tlength\n")
    for record in FastaIterator.iterate(inf):
        scaffold_length = len(list(record.sequence))
        outf.write("%s\t%i\n" % (record.title, scaffold_length))
    outf.close()

############################
# general assembler class
############################
class Assembler(PairedData):

    '''
    general class for assembly algorithms
    '''

    def __init_(self):

        self.kmer = 0
        self.format = None
        self.read_type = None
        self.exp_cov = None
        self.stats_file = None
        self.compressed = False

##########################
# meta-velvet
##########################
class Metavelvet(Assembler):

    '''
    velvet single genome assembly software
    '''

    def build(self, infile, PARAMS):
        '''
        run velveth and velvetg
        followed by meta-velvetg
        '''
        outdir = P.getTempDir(".")
        format = self.getFormat(infile)
        paired = self.checkPairs(infile)

        if not paired:
            pair = ""
            files = infile
            read_type = "short"
        else:
            pair = "-%s" % paired[0]
            files = " ".join([infile, paired[1]])
            read_type = "shortPaired"

        if format == "fastq.1.gz":
            format = "fastq.gz"
        metavelvet_dir = os.path.join(os.getcwd(), "metavelvet.dir")
        track = self.getTrack(os.path.basename(infile))

        self.stats_file = track + ".stats.txt"

        if paired:
            insert_length = "-ins_length %i" % PARAMS["velvetg_insert_length"]
        else:
            insert_length = ""

        # velveth and velvetg have to be run to build hash tables and initial
        # de bruijn graphs
        statement = '''%%(velveth_executable)s %(outdir)s %%(kmer)i -%(format)s -%(read_type)s %(pair)s %(files)s >> %(metavelvet_dir)s/%(track)s_velveth.log
                      ; checkpoint
                      ; mv %(outdir)s/Log %(metavelvet_dir)s/%(track)s.velveth.log
                      ; %%(velvetg_executable)s %(outdir)s -exp_cov auto %(insert_length)s
                      ; checkpoint
                      ; %%(metavelvet_executable)s %(outdir)s %(insert_length)s
                      ; mv %(outdir)s/Roadmaps %(metavelvet_dir)s/%(track)s.roadmaps
                      ; gzip %(metavelvet_dir)s/%(track)s.roadmaps
                      ; mv %(outdir)s/Sequences %(metavelvet_dir)s/%(track)s.sequences
                      ; gzip %(metavelvet_dir)s/%(track)s.sequences
                      ; mv %(outdir)s/Graph2 %(metavelvet_dir)s/%(track)s.graph2
                      ; gzip %(metavelvet_dir)s/%(track)s.graph2
                      ; cat %(outdir)s/meta-velvetg.contigs.fa | python %%(scriptsdir)s/rename_contigs.py -a metavelvet 
                                                                 --log=%(metavelvet_dir)s/%(track)s.contigs.log
                        >  %(metavelvet_dir)s/%(track)s.contigs.fa
                      ; sed -i 's/in/_in/g' %(outdir)s/meta-velvetg.Graph2-stats.txt
                      ; mv  %(outdir)s/meta-velvetg.Graph2-stats.txt %(metavelvet_dir)s/%(track)s.stats.txt
                      ; rm -rf %(outdir)s
                      ''' % locals()
        return statement


##########################
# meta-idba
##########################
class Idba(Metavelvet):

    '''
    meta-idba contig assembler
    '''

    def preprocess(self, infile):
        '''
        fastq files need to be converted to fasta
        and pairs need to be merged
        '''

        mtype = None

        # check for paired end data either in the same file or in a separate file
        # for each read - will need to be gunzipped
        # check compression status
        if infile.endswith(".gz"):
            # check for paired data in separate files
            if len(self.checkPairs(infile)) > 1:
                read1 = infile
                read2 = self.checkPairs(infile)[1]
                temp = P.getTempDir()
            elif self.checkPairs == "interleaved":
                infile_new = os.path.join(temp, P.snip(infile, ".gz"))
                zippy = """gunzip -c %(infile)s > %(infile_new)s; """ % locals()
        else:
            zippy = ""

        # only need to convert if the data are in fastq format
        # reads are fastq and paired in separate files
        if self.getFormat(infile).find("fastq") != -1 and len(self.checkPairs(infile)) > 1:
            mtype = "--merge"  # argument for conversion tool
        # reads are fastq and in the same file
        elif self.getFormat(infile).find("fastq") != -1 and self.checkPairs(infile) == "interleaved":
            mtype = "--paired"  # argument for conversion tool

        # requires a merge of the fastq files in to fasta format
        if mtype:  # the reads are paired end
            if mtype == "--merge":
                outf = P.snip(os.path.basename(read1), ".fastq.1.gz") + ".fa"

                # check if file exists - metaphlan also performs this
                # preprocessing step
                if not os.path.exists(outf):
                    statement = '''python %%(scriptsdir)s/fastqs2fasta.py -a %(read1)s -b %(read2)s --log=%(read1)s.log > %(outf)s
                                ''' % locals()
                    P.run()
                else:
                    E.info("no need to create file %s - exists" % outf)

            elif mtype == "--paired":
                outf = P.snip(os.path.basename(infile_new), ".fastq") + ".fa"
                statement = '''%(zippy)s'''
                P.run()
                statement = '''fq2fa %(mtype)s %(infile_new)s %(outf)s
                             rm -rf %(temp)s''' % locals()
                P.run()
        else:
            statement = None
        return statement

    def build(self, infile):
        '''
        build statement for running idba
        input is merged fasta file and 
        output is contigs fasta file
        '''
        track = self.getTrack(infile)
        outdir = "idba.dir"

        # get temporary file for running idba
        tempdir = P.getTempDir()

        # NB at the moment we assume the default maxkmer of 100
        statement = '''%%(idba_executable)s -r %(infile)s -o %(tempdir)s %%(idba_options)s
                       ; mv %(tempdir)s/scaffold.fa idba.dir/%(track)s.scaffolds.fa
                       ; cat %(tempdir)s/contig-%%(idba_maxkmer)s.fa | python %%(scriptsdir)s/rename_contigs.py -a idba --log=%(outdir)s/%(track)s.contigs.log
                        > idba.dir/%(track)s.contigs.fa''' % locals()

        shutil.rmtree(tempdir)
        return statement


##########################
# Ray meta
##########################
class Ray(Idba):

    '''
    ray contig assembler
    '''

    def build(self, infile):
        '''
        build statement for running Ray
        '''
        track = self.getTrack(os.path.basename(infile))

        format = self.getFormat(infile)
        paired = self.checkPairs(infile)

        tempdir = P.getTempDir(dir=".")

        # check whether the data are paired-end
        if not paired:
            pair = paired
            files = os.path.join(
                tempdir, P.snip(os.path.basename(infile), ".gz"))
            gunzy = "gunzip -c %(infile)s > %(files)s" % locals()
        else:
            pair = paired[0]
            # Ray doesn't like .fastq.1.gz etc
            read1 = infile
            read2 = paired[1]
            read1_new = os.path.join(
                tempdir, read1.replace(".fastq.1.gz", ".1.fastq"))
            read2_new = os.path.join(
                tempdir, read2.replace(".fastq.2.gz", ".2.fastq"))
            files = " ".join([read1_new, read2_new])
            gunzy = """gunzip -c %(read1)s > %(read1_new)s
                       ; gunzip -c %(read2)s > %(read2_new)s""" % locals()

        # ray likes an output directory but needs it not
        # to exist beforehand
        raydir = os.path.join(os.getcwd(), "ray.dir/export_%s" % track)
        raydir_orig = os.path.join(os.getcwd(), "ray.dir")

        # Ray picks up file types so should just have to
        # say whether its paired or not
        # build statement
        common_options = "-k %(kmer)s"
        if pair == "interleaved":
            filetype = "-i"
        elif not pair:
            filetype = "-s"
        elif pair == "separate":
            filetype = "-p"
        else:
            raise IOError, "do not support file of this type: %s" % infile

        # note restrict use to 10 cores

        statement = ''' %(gunzy)s
                       ; mpiexec %%(ray_executable)s %(common_options)s %(filetype)s %(files)s -o %(raydir)s >> %(raydir_orig)s/%(track)s.log
                       ; checkpoint; mv %(raydir)s/Scaffolds.fasta %(raydir_orig)s/%(track)s.scaffolds.fa
                       ; mv %(raydir)s/ScaffoldComponents.txt %(raydir_orig)s/%(track)s.scaffold_components.txt
                       ; mv %(raydir)s/ScaffoldLengths.txt %(raydir_orig)s/%(track)s.scaffold_lengths.txt
                       ; mv %(raydir)s/ScaffoldLinks.txt %(raydir_orig)s/%(track)s.scaffold_links.txt
                       ; cat %(raydir)s/Contigs.fasta | python %%(scriptsdir)s/rename_contigs.py -a ray --log=%(raydir_orig)s/%(track)s.contigs.log
                         > %(raydir_orig)s/%(track)s.contigs.fa
                       ; mv %(raydir)s/OutputNumbers.txt %(raydir_orig)s/%(track)s.numbers.txt
                       ; mv %(raydir)s/CoverageDistribution.txt %(raydir_orig)s/graph/%(track)s.coverage_distribution.txt
                       ; mkdir %(raydir)s/graph
                       ; mv %(raydir)s/CoverageDistributionAnalysis.txt %(raydir_orig)s/graph/%(track)s.coverage_distribution_analysis.txt
                       ; mv %(raydir)s/degreeDistribution.txt %(raydir_orig)s/graph/%(track)s.degree_distribution.txt
                       ; mv %(raydir)s/Kmers.txt %(raydir_orig)s/graph/%(track)s.kmers.txt
                       ; mkdir %(raydir)s/assembly
                       ; mv %(raydir)s/SeedLengthDistribution.txt %(raydir_orig)s/assembly/%(track)s.seed_length_distribution.txt
                       ; mv %(raydir)s/LibraryStatistics.txt %(raydir_orig)s/%(track)s.library_statistics.txt
                       ; mv %(raydir)s/LibraryData.xml %(raydir_orig)s/%(track)s.library_data.xml 
                       ; rm -rf %(tempdir)s''' % locals()

        return statement

##########################
# SGA
##########################
class SGA(Idba):

    '''
    String Graph assembler
    '''

    def build(self, infile):
        '''
        build statement for running SGA
        '''
        track = self.getTrack(os.path.basename(infile))

        # decide which algorithm to use based on
        # read length
        if "%(sga_long)s":
            index_algorithm = "sais"
        else:
            index_algorithm = "ropebwt"

        format = self.getFormat(infile)
        paired = self.checkPairs(infile)

        # directory in which to do the assembly
        tempdir = P.getTempDir(dir=".")

        # check whether the data are paired-end
        if not paired:
            pe_mode = "--pe-mode=0"
            files = os.path.abspath(infile)
        else:
            # DOESN'T DEAL WITH INTERLEAVED FILES YET
            pe_mode = "--pe-mode=1"
            files = " ".join([os.path.abspath(infile), paired[1]])

        executable = "%(sga_executable)s"

        outdir = os.path.abspath("sga.dir")
        ###############################################
        # preprocessing step converts missing bases to
        # random bases or removes sequences with
        # missing bases
        ###############################################
        preprocess_options = "%(sga_preprocess_options)s"
        # outputs a merged fastq file
        outf_preprocessed = track + ".fastq"
        preprocess_statement = "cd %(tempdir)s; %(executable)s preprocess %(pe_mode)s %(preprocess_options)s %(files)s \
                                 -o %(outf_preprocessed)s 2> %(outdir)s/%(track)s_preprocess.log"

        ###############################################
        # indexing reads with FM index
        ###############################################
        index_options = "%(sga_index_options)s"
        index_statement = "%(executable)s index --algorithm=%(index_algorithm)s \
                           %(outf_preprocessed)s 2> %(outdir)s/%(track)s_index.log"

        ###############################################
        # correct sequencing errors in reads
        ###############################################
        correction_method = "%(sga_correction_method)s"

        # if correction_method == "kmer":
        # ADD WARNING HERE
        correction_options = "%(sga_kmer_correction_options)s"
        # elif correction_method == "hybrid":
        #     correction_options = "%(sga_hybrid_correction_options)s"
        # elif correction_method == "overlap":
        #     correction_options = "%(sga_overlap_correction_options)s"
        # else:
        #     raise ValueError("method %s does not exist: choose one of kmer, hybrid, overlap" % correction_method)
        outf_corrected = track + "_corrected.fa"
        metrics = "--metrics=%(track)s.metrics" % locals()
        correction_prefix = os.path.join(
            tempdir, P.snip(outf_corrected, ".fa"))
        correction_statement = "%(executable)s correct %(metrics)s  \
                                --algorithm=%(correction_method)s \
                                %(correction_options)s \
                                %(outf_preprocessed)s \
                                -o %(outf_corrected)s 2> %(outdir)s/%(track)s_corrected.log"

        ###############################################
        # filter low quality reads and low abundance
        # kmers
        ###############################################
        filter_options = "%(sga_filter_options)s"
        outf_filtered = track + "_filtered.fa"
        filter_statement = "sga index %(outf_corrected)s; \
                            %(executable)s filter %(filter_options)s  \
                            -o %(outf_filtered)s \
                            %(outf_corrected)s 2> %(outdir)s/%(track)s_filtered.log"

        ###############################################
        # overlap reads
        ###############################################
        # Note "asqg" is the default output from sga
        outf_overlap = track + "_filtered.asqg.gz"
        threads = "%(sga_threads)s"
        overlap_options = "%(sga_overlap_options)s"
        overlap_statement = "%(executable)s overlap %(overlap_options)s \
                             %(outf_filtered)s 2> %(outdir)s/%(track)s_overlap.log"

        ###############################################
        # assemble reads and perform error removal
        ###############################################
        assembly_options = "%(sga_assembly_options)s"
        error_removal_options = "%(sga_error_removal_options)s"
        out_prefix = track
        assembly_statement = "%(executable)s assemble %(assembly_options)s \
                              %(outf_overlap)s \
                              --out-prefix=%(out_prefix)s 2> %(outdir)s/%(track)s_contigs.log"

        ###############################################
        # build statement
        ###############################################
        metrics_file = os.path.basename(metrics.replace("--metrics=", ""))
        contigs_file = os.path.basename(out_prefix + "-contigs.fa")
        move_statement = "mv %(metrics_file)s %(outdir)s/%(metrics_file)s; \
                          cat %(contigs_file)s \
                          | python %%(scriptsdir)s/rename_contigs.py \
                          --log=%(outdir)s/%(track)s.contigs.log \
                          -a sga > %(outdir)s/%(track)s.contigs.fa"

        statement = "; ".join([preprocess_statement,
                               index_statement,
                               correction_statement,
                               filter_statement,
                               overlap_statement,
                               assembly_statement,
                               move_statement,
                               "rm -rf %(tempdir)s"]) % locals()

        return statement

##########################
# SOAPdenovo2
##########################
class SoapDenovo2(Idba):

    '''
    soapdenovo is a single genome assembler
    and therefore may belong elsewhere
    '''

    def config(self, infile, outfile, PARAMS):
        '''
        build the configuration file based on
        parameters specified in the .ini file
        '''
        # output directory
        outdir = "soapdenovo.dir"

        # config file
        outf = open(outfile, "w")
        paired = self.checkPairs(infile)

        # global parameters - only specifies maximum
        # read length at the moment - shorter are
        # trimmed
        max_rd_len = "%(soapdenovo_max_rd_len)s" % PARAMS

        # assume library parameters are the same
        # for each library
        avg_ins = "%(soapdenovo_avg_ins)s" % PARAMS
        asm_flags = "%(soapdenovo_asm_flags)s" % PARAMS
        reverse_seq = "%(soapdenovo_reverse_seq)s" % PARAMS
        rank = "%(soapdenovo_rank)s" % PARAMS

        # specify files
        if not paired:
            q = "q=" + infile
        else:
            q = "\n".join(["q1=" + infile, "q2=" + paired[1]])
        outf.write(
            """max_rd_len=%(max_rd_len)s\n[LIB]\navg_ins=%(avg_ins)s\nasm_flags=%(asm_flags)s\nreverse_seq=%(reverse_seq)s\nrank=%(rank)s\n%(q)s\n""" % locals())

        outf.close()

    def build(self, config):
        '''
        return build statement to be run
        '''
        # output directory
        outdir = "soapdenovo.dir"

        # get track from config file
        for line in open(config).readlines():
            if line.startswith("q") or line.startswith("q1"):
                track = self.getTrack(line[:-1].split("=")[1])

        options = "%(soapdenovo_options)s"
        tempdir = P.getTempDir(".")
        statement = '''%%(soapdenovo_executable)s all 
                       -s %%(infile)s 
                       -o %(tempdir)s/%(track)s
                       -K %%(kmer)s
                       %(options)s; checkpoint;
                       mv %(tempdir)s/%(track)s* %(outdir)s;
                       mv %(outdir)s/%(track)s.contig %(outdir)s/%(track)s.contigs.fa;
                       rm -rf %(tempdir)s''' % locals()
        return statement


##########################
# Spades
##########################
class Spades(Idba):

    def build(self, infile):
        '''
        build statement for running spades
        '''
        track = self.getTrack(os.path.basename(infile))

        format = self.getFormat(infile)
        paired = self.checkPairs(infile)

        tempdir = P.getTempDir(".")
        outdir = "spades.dir"

        # input files
        if not paired:
            files = infile
            files_statement = "-s %s" % files
        else:
            files_statement = "-1 "+ "-2 ".join( [infile, paired[1]] )

        # kmer to use
        k = "-k %(kmer)s"
        
        # spades options
        spades_options = "%(spades_options)s"

        # deal with spades output
        move_statement = """mv %(tempdir)s/corrected/%(track)s*cor.* %(outdir)s; \
                          mv %(tempdir)s/contigs.fasta %(outdir)s/%(track)s.contigs.fa; \
                          mv %(tempdir)s/scaffolds.fasta %(outdir)s/%(track)s.scaffolds.fa \
                          mv %(tempdir)s/spades.log %(outdir)s/%(track)s.contigs.log""" % locals()
        
        # statement - simple and default
        statement = '''spades.py %(files_statement)s -o %(tempdir)s %(k)s %(spades_options)s;
                       %(move_statement)s;
                       rm -rf %(tempdir)s''' % locals()
        return statement


##########################
# metaphlan
##########################
class Metaphlan(Idba):

    '''
    metphlan is more of an annotation tool
    and therefore may be removed from this pipeline
    - however it is directly relevant for metagenome sequencing
    '''

    def build(self, infile, method="read_map"):
        '''
        build statement for running metaphlan
        '''
        track = self.getTrack(infile)

        if method == "read_map":
            statement = '''cat %(infile)s                                                                                                                                                                                                           
                      | metaphlan.py -t reads_map                                                                                                                                                                              
                       --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map                                                                                                                                                                      
                      | python %%(scriptsdir)s/metaphlan2table.py -t read_map                                                                                                                                                              
                       --log=%%(outfile)s.log                                                                                                                                                                                                     
                      > %%(outfile)s; checkpoint
                      ; sed -i 's/order/_order/g' %%(outfile)s''' % locals()
        elif method == "rel_ab":
            statement = '''cat %(infile)s                                                                                                                                                                                                           
                      | metaphlan.py -t rel_ab                                                                                                                                                                              
                       --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map                                                                                                                                                                                                                                                                                                                                   
                      | python %%(scriptsdir)s/metaphlan2table.py -t rel_ab 
                       --log=%%(outfile)s.log                                                                                                                                                                                                   
                      > %%(outfile)s''' % locals()
        else:
            raise ValueError, "do not support method %s" % method

        return statement

         # DEPRECATED
        # create single fasta file if required (reads are fastq format)
        # inf = track + ".fa"

        # if not os.path.exists(inf):
        #     job_options = " -l mem_free=30G"
        #     to_cluster = True
        #     if self.preprocess(infile):
        #         statement = self.preprocess(infile)
        #         P.run()
        #     else:
        #         inf = infile

        # if method == "read_map":
        #     statement = '''cat %(inf)s
        #               | python %%(scriptsdir)s/metaphlan.py -t reads_map
        #                --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map
        #               | python %%(scriptsdir)s/metaphlan2table.py -t read_map
        #                --log=%%(outfile)s.log
        #               > %%(outfile)s; checkpoint
        #               ; sed -i 's/order/_order/g' %%(outfile)s''' % locals()
        # elif method == "rel_ab":
        #     statement = '''cat %(inf)s
        #               | python %%(scriptsdir)s/metaphlan.py -t rel_ab
        #                --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map
        #               | python %%(scriptsdir)s/metaphlan2table.py -t rel_ab
        #                --log=%%(outfile)s.log
        #               > %%(outfile)s''' % locals()
        # else:
        #     raise ValueError, "do not support method %s" % method

        # return statement


##########################
# cortex_var
##########################
class Cortex_var(Idba):

    '''
    cortex genome assembler
    '''

    def build(self, infile):

        track = self.getTrack(infile)

        format = self.getFormat(infile)
        if format.endswith(".gz"):
            format = P.snip(format, ".gz")
        format = format.upper()

        # cortex_var only uses paired end information to
        # remove pcr duplicates
        if not self.checkPairs(infile):
            paired = "--se_list"
            reads = os.path.join(os.getcwd(), infile)

        elif len(self.checkPairs(infile)) > 1:
            paired = "--pe_list"
            read1 = infile
            format = P.snip(format, ".1")
            read2 = self.checkPairs(infile)[1]

        elif self.checkPairs(infile) == "interleaved":
            raise ValueError, "pipeline does not support file of type 'interleaved'"

        temp = P.getTempDir()
        read1_new = os.path.join(temp, P.snip(read1, ".1.gz"))
        read2_new = os.path.join(temp, P.snip(read2, ".2.gz"))

        # paired end list
        list1 = open("cortex_var.dir/read1.txt", "w")
        list2 = open("cortex_var.dir/read2.txt", "w")
        list1.write(read1_new + "\n")
        list2.write(read2_new + "\n")
        list1.close()
        list2.close()

        list1 = os.path.abspath("cortex_var.dir/read1.txt")
        list2 = os.path.abspath("cortex_var.dir/read2.txt")

        reads = ",".join([os.path.join(os.getcwd(), x)
                         for x in [read1_new, read2_new]])
        statement = '''  gunzip -c %(read1)s > %(read1_new)s
                       ; gunzip -c %(read2)s > %(read2_new)s  
                       ; cd cortex_var.dir
                       ; %%(cortex_var_executable)s %(paired)s %(list1)s,%(list2)s 
                       --format %(format)s
                       --mem_height 15
                       --quality_score_threshold %%(cortex_var_qual_threshold)i 
                       --remove_pcr_duplicates 
                       --remove_low_coverage_supernodes %%(cortex_var_rm_low_coverage_supernodes)i
                       --sample_id %(track)s
                       --kmer_size %%(kmer)s
                       --dump_binary dump_binary.ctx
                       ; rm -rf %(temp)s
                    ''' % locals()

        return statement
