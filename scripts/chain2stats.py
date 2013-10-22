'''
chain2stats.py
==============

:Author: Stephen Sansom
:Release: $Id$
:Date: |today|
:Tags: Genomics GenomeAlignment Summary CHAIN

Purpose
-------

This script computes the coverage of the target and query genomes
represented by a UCSC chain file.

Usage
-----

Example::

   python chain2stats.py < sacCer3ToSacCer1.over.chain.gz

Type::

   python chain2stats.py --help

for command line help.

Documentation
-------------

For large chain file this script requires a substantial amount of
memory for chromosome pair analysis.

Command line options
---------------------
'''

#import modules
import os
import sys
from operator import add
from numpy import *
from operator import itemgetter,attrgetter
import bx.bitset
import bx.bitset_builders
import collections

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

############################ Functions/Generators ############################

def chain_iterator(infile):
    lines = []
    for line in infile:
        if line.startswith("#"): continue
        if line.strip() == "": continue 
        if line.startswith("chain"):
            if lines: yield lines
            lines = []
        lines.append( line )
    yield lines

######################### Functions/Generators end ############################


###############################################################################
################################ Classes ######################################
###############################################################################

##-----------------------------------------------------------------------------

class Chain:

    "class with methods for getting gapped and ungapped regions represented by a chain"

    def __init__(self,data):
        self.lines = data
        self.header = self.lines[0]
        self.atts = self.header.split()
        a = len(self.atts)
        if a != 13: #sanity check, can remove
            raise Exception("Unexpected Header length\n"+self.header)
        self.score = int(self.atts[1])
        self.tname = self.atts[2]
        self.tsize = int(self.atts[3])
        self.tstrand = self.atts[4]
        self.tstart = int(self.atts[5])
        self.tend = int(self.atts[6])
        self.qname = self.atts[7]
        self.qsize = int(self.atts[8])
        self.qstrand = self.atts[9]
        b = self.qstrand
        if (b != "+") and (b != "-"):
            raise Exception("unexpected strand id")
        self.qstart = int(self.atts[10])
        self.qend = int(self.atts[11])
        self.id = int(self.atts[12])
        self.tlength = self.tend - self.tstart
        self.qlength = self.qend - self.qstart
        self.nlines = len(self.lines)
        self.nblocks = self.nlines - 1
        (self.tgr,self.qgr) = self._gapped_regions()
        (self.tugr,self.qugr) = self._ungapped_regions()

    def _gapped_regions(self):
        tr,qr = [],[]
        tr.append([self.tstart,self.tlength])
        qr.append([self.qstart,self.qlength])
        return(tr,qr)

    def _ungapped_regions(self):
        ts = self.tstart
        qs = self.qstart
        tr,qr = [],[]
        for i in range(1,self.nlines):
            dataline = self.lines[i]
            values = dataline.split()
            len_align = int(values[0])
            te = ts + len_align
            qe = qs + len_align
            tr.append([ts,len_align])
            qr.append([qs,len_align])
            if (i < self.nlines-1):
                ts = te + int(values[1])
                qs = qe + int(values [2])
        return(tr,qr)

    def get_regions(self,gapped):
        if gapped == True:
            treg, qreg = self.tgr, self.qgr
        else:
            treg, qreg = self.tugr, self.qugr            
        return(treg,qreg)

##-----------------------------------------------------------------------------

class Bitsets_Container:

    "General class for holding bitsets in a {key1: {key2: bitset }} structure"

    def __init__(self):
        self.bitsets = collections.defaultdict( dict )

    def addRanges(self,key1,key2,size,ranges):
        if key2 not in self.bitsets[key1]:
            if size == 199501827: # this randomly segfaults a binnedbitset
                self.bitsets[key1][key2] = bx.bitset.BitSet(size)
            else:
                self.bitsets[key1][key2] = bx.bitset.BinnedBitSet(size) 
        for r in ranges: #bx-python does sanity checks.
           self.bitsets[key1][key2].set_range(r[0],r[1]) 

    def bitset_coverage(self,key1,key2):
        b = self.bitsets[key1][key2]
        l = b.size
        c = b.count_range( 0, l )
        return (c , l)

##**---------------------------------------------------------------------------

class ChainCounter():

    '''
    Top, defining class for a chain counter.
    
    Counters must have:
    (1) an add method to add a new chain
    (2) a write_report method

    Counters will generally initialise a container and add metrics to it when passed a chain
    '''

#---------- abstract methods------------

    def add(self,c):
        raise NotImplementError("abstract method - implement in subclass")

    def _get_stats(self,options):
        '''Calculates and stores the stats if necessary - called by all report methods'''
        raise NotImplementError("abstract method - implement in subclass")

    def report(self,options):
        '''should call the "write_report" method to do the writing'''
        raise NotImplementError("abstract method - implement in subclass")

    def tabbed_report(self,options):
        '''should call the "write_tabbed" method to do the writing'''
        raise NotImplementError("abstract method - implement in subclass")


#----------- methods --------------------

    def _get_basic_stats(self,numlist,string="{:.2f}"):
        s = [string.format(mean(numlist)),string.format(median(numlist)),string.format(max(numlist)),string.format(min(numlist))]
        return(s)

    def _wrap_basic_stats(self,s,tabbed=False):
        if tabbed==False:
            r = ''.join(["Mean:",s[0],", median:",s[1],", max:",s[2],", min:",s[3],"\n"])
        else: 
            r = '\t'.join(s)
        return(r)

    def _write_report(self,options,report):
        options.stdout.write(''.join(['\n','\n'.join(report),'\n']))
        
    def _write_tabbed(self,name,lines,E):
        outfile = E.openOutputFile(name)
        outfile.write('\n'.join(lines))
        outfile.write('\n')
        outfile.close

    def _wrap_header(self,text=""):
        '''decorates header for report'''
        if text == "": text = self.header
        underline = ''.join(["-" for i in text])
        return([text,underline])

##----------------------------------------------------------------------------

class CounterPerChromosome(ChainCounter):

    header = "Report for chains per chromosome"
    
    def __init__(self,gapped):
        self.container = Bitsets_Container()
        self.gapped = gapped
        self.cov = 0
        if self.gapped == True: self.name="gapped"
        else: self.name = "ungapped"

    def _pc_coverage(self,coverage,total_length):
        pc_cov = float(coverage)/float(total_length) * 100
        return(pc_cov)

    def _coverage(self,key1):
        total, total_len = 0,0
        d = self.container.bitsets[key1]
        for key2 in sorted(d):
            (s,l) = self.container.bitset_coverage(key1,key2)
            total += s
            total_len += l
        gen_coverage = self._pc_coverage(total,total_len)
        return gen_coverage

    def add(self,c):
        treg, qreg = c.get_regions(self.gapped)
        self.container.addRanges("target",c.tname,c.tsize,treg)
        self.container.addRanges("query",c.qname,c.qsize,qreg)

    def _get_stats(self):
        if self.cov == 0: 
            self.cov = {}
            for i in ("target","query"):
                self.cov[i] = self._coverage(i)

    def report(self,options):
        self._get_stats()
        report = self._wrap_header("".join([self.header," ",self.name," statistics"]))
        report.append("Target genome {:.2f}, query genome: {:.2f}\n".format(self.cov["target"],self.cov["query"]))
        self._write_report(options,report)

    def tabbed_report(self,options,E):
        self._get_stats()
        lines = ["target_genome\tquery_genome"]
        lines.append('{:.2f}\t{:.2f}'.format(self.cov["target"],self.cov["query"]))
        self._write_tabbed(''.join(["genome_coverage_",self.name]),lines,E)        

##-----------------------------------------------------------------------------

class CounterPerChromosomePair(CounterPerChromosome):
    
    header = "Report for chains per chromosome pair"

    def _make_table(self,options):
        table = []
        table.append("t_chrom\tt_cov\tq_chrom\tq_cov\tbases aligned")
        tchrom = self.cov[0][0]
        j = 0
        for i in self.cov:
            if i[0]!=tchrom:
                tchrom = i[0]
                j=0       
            if j < options.nperchrom:
                table.append("{}\t{:.2f}\t{}\t{:.2f}\t{:,}".format(i[0],i[1],i[2],i[3],i[4]))
                j += 1
            else:
                continue
        table.append("")
        return(table)

    def add(self,c):
        treg, qreg = c.get_regions(self.gapped)
        if c.qstrand == "+": strand = "plus" 
        else: strand = "minus"
        pair_name = (c.tname, c.qname)
        self.container.addRanges(pair_name,"target",c.tsize,treg)
        self.container.addRanges(pair_name,"query",c.qsize,qreg)
        self.container.addRanges(pair_name,strand,c.qsize,qreg)

    def _get_stats(self,options):
        if self.cov == 0:
            result = []
            for key1 in self.container.bitsets:
                (tm,tl) = self.container.bitset_coverage(key1,'target')
                tc = self._pc_coverage(tm,tl)
                (qm,ql) = self.container.bitset_coverage(key1,'query')
                qc = self._pc_coverage(qm,ql)
                result.append([key1[0],tc,key1[1],qc,tm])
            a = sorted(result,key=itemgetter(4),reverse=True)
            self.cov = sorted(a,key=itemgetter(0))
            self.cov_table = self._make_table(options)

    def report(self,options):
        self._get_stats(options)
        report = (self._wrap_header("".join([self.header," ",self.name," statistics"])))
        report = report + self.cov_table
        self._write_report(options,report)

    def tabbed_report(self,options,E):
        self._get_stats(options)
        self._write_tabbed(''.join(["chrom_pair_cov_",self.name]),self.cov_table,E)        

##-----------------------------------------------------------------------------

class CounterOfGappedChainLengths( ChainCounter ):
    
    header = "Report for chain lengths"

    def __init__(self,gapped):
        self.nchains = 0
        self.tcls,self.qcls = [],[]
        self.stats = {}
        self.gapped = gapped
        if self.gapped == True: self.name="gapped"
        else: self.name = "ungapped"

    def _chain_stats_report(self,name):
        report = [''.join(["Stats for ",name," chain lengths:"])]
        report.append(self._wrap_basic_stats(self.stats[name]))
        return(report)      

    def add(self,c):
        treg, qreg = c.get_regions(self.gapped)
        n = len(treg)
        if n != len(qreg):
            raise Exception("Different numbers of chains for query and target!?!")
        for i in range(0,n):
            self.tcls.append(treg[i][1])
            self.qcls.append(qreg[i][1])
            self.nchains += 1

    def _get_stats(self):
        if ("target" in self.stats)==0:
            self.stats["target"] = self._get_basic_stats(self.tcls,string="{:,.0f}")
        if ("query" in self.stats)==0:
            self.stats["query"] = self._get_basic_stats(self.qcls, string="{:,.0f}")

    def report(self,options):
        self._get_stats()
        report = (self._wrap_header(" ".join([self.header,self.name])))
        report.append("number of chains = {:}\n".format(self.nchains))
        for i in ("target","query"): 
            report = report + self._chain_stats_report(i)
        self._write_report(options,report)

    def tabbed_report(self,options,E):
        self._get_stats()
        for i in ("target","query"):
            lines = ["mean\tmedian\tmax\tmin"]
            lines.append(self._wrap_basic_stats(self.stats[i],tabbed=True))
            self._write_tabbed('_'.join([i,"lengths",self.name]),lines,E)        

##-----------------------------------------------------------------------------

class CounterPercentIdentify( ChainCounter ):

    header = "Report on Percent Indentities"

    def __init__(self,tname,qname):
        self.tfasta = IndexedFasta(tname)
        self.qfasta = IndexedFasta(qname)
        self.pids = []
        self.stats = 0

    def _get_pid(self,x,y):
        z = zip(x,y)
        pid = (float(len([x for x,y in z if x==y])) / float(len(z)) * 100)
        return(pid)

    def add(self,c):
        nreg = len(c.tugr)
        for i in range(0,nreg):                                                            
            tseq = self.tfasta.getSequence(c.tname,"+",c.tugr[i][0],(sum(c.tugr[i])))
            qseq = self.qfasta.getSequence(c.qname,c.qstrand,c.qugr[i][0],(sum(c.qugr[i])))
            pid = self._get_pid(tseq.lower(),qseq.lower())
        self.pids.append( pid )
        
    def _get_stats(self):
        if self.stats == 0:
            self.stats = self._get_basic_stats(self.pids,string="{:.2f}")            

    def report(self,options):
        self._get_stats()
        report = self._wrap_header()
        report.append(self._wrap_basic_stats(self.stats))
        self._write_report(options,report)

    def tabbed_report(self,options,E):
        self._get_stats()
        lines = ["mean\tmedian\tmax\tmin"]
        lines.append(self._wrap_basic_stats(self.stats,tabbed=True))
        self._write_tabbed("pids",lines,E)        
       
##-----------------------------------------------------------------------------

class CounterOfErrors( ChainCounter ):

    '''class for reporting invalid contig sizes in chains'''

    header = "Contig size validation report"

    def  __init__(self,options):
        self.tdb = IndexedFasta(options.dbpath+options.targetgenome)
        self.qdb = IndexedFasta(options.dbpath+options.querygenome)
        self.tcontigs = self.tdb.getContigSizes()
        self.qcontigs = self.qdb.getContigSizes()
        self.badchains = []

    def add(self,c):
        db_tsize = self.tcontigs[c.tname]
        db_qsize = self.qcontigs[c.qname]
        if c.tsize != db_tsize:
            self.badchains.append('\t'.join([str(x) for x in c.atts] + [" #bad target contigsize"]))
        if c.qsize != db_qsize:
            self.badchains.append('\t'.join([str(x) for x in c.atts] + [" #bad query contigsize"]))
    
    def report(self,options):
        report = self._wrap_header()
        if len(self.badchains) == 0:
            report.append("All chains passed validation")
        else:
            report = report + self.badchains
        self._write_report(options,report)

    def tabbed_report(self,options,E):
        if len(self.badchains) > 0:
            lines = self.badchains
        else:
            lines = ["#no bad chains found"]
        self._write_tabbed("bad_contig_sizes",lines,E)        

##-----------------------------------------------------------------------------

###############################################################################
############################# Classes end #####################################
###############################################################################


def main(argv = None ):
    if not argv: argv = sys.argv
    
    #get the options
    parser = E.OptionParser(version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                            usage = globals()["__doc__"] )

    parser.add_option("-c", "--chainfile", dest="chainfile", type="string",
                      help="the chain file to analyse", metavar="FILE")
    parser.add_option("-p", "--perchrom", dest="perchrom", action="store_true",
                      help="Set to 1 or \"True\" to perform per chromosome pair analysis",default=False)
    parser.add_option("-n", "--nperchrom", dest="nperchrom", type="int",
                      help="Number of aligments to report on per chromosome pair",default=2)
    parser.add_option("-i","--identity", dest="identity", action="store_true",
                      help="Generate stats on the sequence identity of the gapped chains. Requires FastIndex.py", default=False)
    parser.add_option("-d", "--dbpath", dest="dbpath", type="string",
                      help="The path to the indexed fasta files",default="/ifs/mirror/genomes/plain/")
    parser.add_option("-t", "--targetgenome", dest="targetgenome", type="string",
                      help="The target genome, eg. Mm19",default=False)
    parser.add_option("-q", "--querygenome", dest="querygenome", type="string",
                      help="The query genome eg. Hg17",default=False)
    parser.add_option("-e", "--errors", dest="errors", action="store_true",
                      help="Check chains for erroneous contig sizes using the given db",default=False)
    parser.add_option("-r", "--report", dest="report", action="store_true",
                      help="Write out tab-delimited reports for each analysis",default=False)

    (options, args) = E.Start(parser, argv = argv, add_output_options = True )
    
    #make a list of counting objects
    counters = []

    counters.append(CounterPerChromosome( gapped = True ) )
    counters.append(CounterPerChromosome( gapped = False ) )

    if options.perchrom:
        counters.append( CounterPerChromosomePair( gapped = True ) )
        counters.append( CounterPerChromosomePair( gapped = False ) )

    counters.append(CounterOfGappedChainLengths( gapped = True ))
    counters.append(CounterOfGappedChainLengths( gapped = False ))

    if options.identity == True:
        if options.targetgenome == 0 or options.querygenome == 0:
            raise Exception("Target and query database must be specified with the \"-e\" flag") 
        t_db_path = os.path.join(options.dbpath,options.targetgenome)
        q_db_path = os.path.join(options.dbpath,options.querygenome)
        counters.append(CounterPercentIdentify( t_db_path, q_db_path ))

    if options.errors == True:
        if options.targetgenome == 0 or options.querygenome == 0:
            raise Exception("Target and query database must be specified with the \"-e\" flag") 
        counters.append(CounterOfErrors(options))

    #iterate over the chains and counters
    for chain in chain_iterator(options.stdin):
        c = Chain(chain)
        for counter in counters:
            counter.add(c)
    
    #write a report to stdout and individual reports to tab delimited files
    options.stdout.write("\n\n********** chain2stats report starts **********\n")

    for counter in counters:
        counter.report(options)
        if options.report == True:
            counter.tabbed_report(options,E)

    options.stdout.write("\n********** chain2stats report ends **********\n\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
