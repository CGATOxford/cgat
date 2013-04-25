'''
classes and utility functions for pipeline_genomeassembly.py

Different assembly tools will use different inputs. Some can take
fasta files whereas others will take fastq and in either case can 
be paired-end (in the same or different files) or single end
'''

import sys, re, os, tempfile, collections, shutil, gzip, sqlite3

import IOTools
import Pipeline as P
import Experiment as E
import PipelineMapping
import FastaIterator
import Fastq
import glob
import collections
import PipelineTracks
import metaphlan_utils
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
        possible_formats = ["fasta", "fastq", "fasta.gz", "fastq.gz", "fasta.1.gz", "fastq.1.gz"]
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
        return P.snip(infile, ".%s" % self.fileFormat(infile))

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
        assert len(glob.glob(read2)) > 0, "cannot find %s file for read 2 in the pair for %s" % (read2, infile)
        return ("separate", read2)

    def checkPairs(self, infile):
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
        
        for record in iterator(inf):
            # make sure there are not other "/" in the sequence name
            seq_id = record.title.split("/")
            assert len(seq_id) == 2, "cannot deal with this sequence name %s" % record.title
            seq_id = seq_id[0]
            if seq_id not in pairs:
                pairs.add(seq_id)
            else:
                print "found pair for %s: %s" % (seq_id, record.title)
                paired = "interleaved"
                break
        return paired

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
    for length in sorted(scaffold_lengths, reverse = True):
        while cum_length <= total_length*(float(N)/100):
            index += 1
            cum_length += length

    # output the results                                                                                                                                                                                                                     
    outf = open(stats_file, "w")
    outf.write("nscaffolds\tscaffold_length\tN%i\tmedian_length\tmean_length\tmax_length\n" % N)
    outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (number_of_scaffolds, total_length, sorted(scaffold_lengths, reverse = True)[index], str(median_length), str(mean_length), str(max_length)))

###############################
###############################
###############################
def build_scaffold_lengths(contigs_file, outfile, params):
    '''
    output the distribution of scaffold lengths
    '''
    PARAMS = params

    if PARAMS["filter"]:
        f = PARAMS["filter"]
    else:
        f = 0
    inf = open(contigs_file)
    outf = open(outfile, "w")
    outf.write("scaffold_name\tlength\n")
    for record in FastaIterator.iterate(inf):
        scaffold_length = len(list(record.sequence))
        if scaffold_length > f:
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
    def build(self, infile):
        '''
        run velveth and velvetg
        followed by meta-velvetg
        '''
        outdir = P.getTempDir()
        format = self.getFormat(infile)
        paired = self.checkPairs(infile)
        if len(paired) > 1:
            pair = paired[0]
            files = " ".join([infile, paired[1]])
        else:
            pair = paired
            files = infile
        if format == "fastq.1.gz":
            format = "fastq.gz"
        metavelvet_dir = os.path.join(os.getcwd(), "metavelvet.dir")
        track = self.getTrack(infile)
        
        self.stats_file = track + ".stats.txt"

        # velveth and velvetg have to be run to build hash tables and initial de bruijn graphs
        statement = '''%%(velveth_executable)s %(outdir)s %%(kmer)i -%(format)s -shortPaired -%(pair)s %(files)s
                      ; cd %(outdir)s; %%(velvetg_executable)s %(outdir)s -exp_cov auto -ins_length %%(velvetg_insert_length)i
                      ; %%(metavelvet_executable)s %(outdir)s -ins_length %%(velvetg_insert_length)i
                      ; mv %(outdir)s/Roadmaps %(metavelvet_dir)s/%(track)s.roadmaps
                      ; gzip %(metavelvet_dir)s/%(track)s.roadmaps
                      ; mv %(outdir)s/Sequences %(metavelvet_dir)s/%(track)s.sequences
                      ; gzip %(metavelvet_dir)s/%(track)s.sequences
                      ; mv %(outdir)s/Graph2 %(metavelvet_dir)s/%(track)s.graph2
                      ; gzip %(metavelvet_dir)s/%(track)s.graph2
                      ; mv %(outdir)s/meta-velvetg.contigs.fa %(metavelvet_dir)s/%(track)s.contigs.fa
                      ; sed -i 's/in/_in/g' %(outdir)s/meta-velvetg.Graph2-stats.txt
                      ; mv  %(outdir)s/meta-velvetg.Graph2-stats.txt %(metavelvet_dir)s/%(track)s.stats.txt''' % locals()
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
            if len(self.checkPairs(infile)) > 1: # check for paired data in separate files
               read1 = infile 
               read2 = self.checkPairs(infile)[1]
               temp = P.getTempDir()
               read1_new = os.path.join(temp, P.snip(infile, ".gz"))
               read2_new = os.path.join(temp, P.snip(self.checkPairs(infile)[1], ".gz"))
               zippy = """gunzip -c %(read1)s > %(read1_new)s
                       ; gunzip -c %(read2)s > %(read2_new)s; """ % locals()
            elif self.checkPairs == "interleaved":
                infile_new = os.path.join(temp, P.snip(infile, ".gz")) 
                zippy = """gunzip -c %(infile)s > %(infile_new)s; """ % locals()
        else:
            zippy = ""
        
        # only need to convert if the data are in fastq format
        if self.getFormat(infile).find("fastq") != -1 and len(self.checkPairs(infile)) >1: # reads are fastq and paired in separate files
            mtype = "--merge" # argument for conversion tool
        elif self.getFormat(infile).find("fastq") != -1 and self.checkPairs(infile) == "interleaved": # reads are fastq and in the same file
            mtype = "--paired" # argument for conversion tool

        # build statement
        if mtype: # the reads are paired end
            if mtype == "--merge":
                outf = P.snip(os.path.basename(read1_new), ".fastq.1") + ".fa" 
                statement = '''%(zippy)s
                             fq2fa %(mtype)s %(read1_new)s %(read2_new)s %(outf)s
                             ''' % locals()
            elif mtype == "--paired":
                outf = P.snip(os.path.basename(infile_new), ".fastq") + ".fa"
                statement = '''%(zippy)s
                             fq2fa %(mtype)s %(infile_new)s %(outf)s
                             ''' % locals()
        else:
            statement = None
        return statement

    def build(self, infile):
        '''build statement for running idba
        '''
        track = self.getTrack(infile)
        
        # create single fasta file if required (reads are fastq format)
        if self.preprocess(infile):
            inf = track + ".fa"
            if not os.path.exists(inf):
                statement = self.preprocess(infile)
                job_options = " -l mem_free=30G"
                to_cluster = True
                P.run()
        else:
            inf = infile
        
        # build statement
        track = self.getTrack(infile)
        outdir = "idba.dir"
        statement = '''%%(idba_executable)s -r %(inf)s -o %(outdir)s
                       ; mv idba.dir/scaffold.fa idba.dir/%(track)s.scaffold.fa''' % locals()
        return statement
        

##########################
# Ray meta
##########################
# class Ray(Idba):
#     '''
#     ray contig assembler
#     '''
#     def build(self, infile):
#         '''
#         build statement for running Ray
#         '''
#         track = self.getTrack(infile)

#         outdir = P.getTempDir()
#         format = self.getFormat(infile)
#         paired = self.checkPairs(infile)

#         # check whether the data are paired-end
#         if len(paired) > 1:
#             pair = paired[0]
#             files = " ".join([infile, paired[1]])
#         else:
#             pair = paired
#             files = infile
            
#         outdir = os.path.join(os.getcwd(), "ray.dir")
        
#         # TODO Ray picks up file types so should just have to
#         # say whether its paired or not
        
#         # build statement
#         common_options = "-k %%(kmer)s "
#         if pair:
#             statement = '''mpiexec -n 80 %%(ray_executable)s -k %%(kmer)s'''
                         
                         
        

##########################
# metaphlan
##########################
class Metaphlan(Idba):
    '''
    metphlan is more of an annotation tool
    and therefore may be removed from this pipline
    - however it is directly relevant for metagenome sequencing
    '''
    def build(self, infile, method="read_map"):
        '''
        build statement for running metaphlan
        '''
        track = self.getTrack(infile)
        
        # create single fasta file if required (reads are fastq format)
        inf = track + ".fa"
        if not os.path.exists(inf):
            job_options = " -l mem_free=30G"
            to_cluster = True
            if self.preprocess(infile):
                statement = self.preprocess(infile)
                P.run()
            else:
                inf = infile
        
        if method == "read_map":
            statement = '''cat %(inf)s                                                                                                                                                                                                           
                      | python %%(scriptsdir)s/metaphlan.py -t reads_map                                                                                                                                                                              
                       --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map                                                                                                                                                                      
                      | python %%(scriptsdir)s/metaphlan2table.py -t read_map                                                                                                                                                              
                       --log=%%(outfile)s.log                                                                                                                                                                                                   
                      > %%(outfile)s''' % locals()
        elif method == "rel_ab":
            statement = '''cat %(inf)s                                                                                                                                                                                                           
                      | python %%(scriptsdir)s/metaphlan.py -t rel_ab                                                                                                                                                                              
                       --input_type multifasta %%(method)s %%(metaphlan_db)s --no_map                                                                                                                                                                                                                                                                                                                                   
                      | python %%(scriptsdir)s/metaphlan2table.py -t rel_ab 
                       --log=%%(outfile)s.log                                                                                                                                                                                                   
                      > %%(outfile)s''' % locals()
        else:
            raise ValueError, "do not support method %s" % method

        return statement

##########################
# cortex_var 
##########################
class Cortex_var(Idba):
    '''
    cortex genome assembler
    '''
    def build(self,infile):
        
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

        reads = ",".join([os.path.join(os.getcwd(), x) for x in [read1_new, read2_new]])        
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
                    ''' % locals()
        
        return statement
    







