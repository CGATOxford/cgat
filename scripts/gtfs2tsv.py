'''
gtfs2tsv.py - compare two genesets
==================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares two genesets in :term:`gtf`-formatted files and output lists
of shared and unique genes.

It outputs the results of the comparison into various sections. The sections are
split into separate output files whose names are determined by the ``--output-filename-pattern`` 
option. The sections are:

genes_ovl
   Table with overlapping genes

genes_total
   Summary statistic of overlapping genes

genes_uniq1
   List of genes unique in set 1

genes_uniq2
   List of genes unique in set 2

Usage
-----

Example::

   python gtfs2tsv.py a.gtf b.gtf

Type::

   python gtfs2tsv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import re
import os
import optparse
import collections
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import bx.intervals.intersection

def GetNextLine( infile ):
    for line in infile:
        if line[0] == "#": continue
        return line
    return None

class Counts:

    mPercentFormat = "%.2f"

    def __init__(self, add_percent ):

        self.nleft, self.nright, self.noverlap = 0, 0, 0
        self.nunique_left = 0
        self.nunique_right = 0
        self.nidentical = 0
        self.nhalf = 0
        self.nsplit_left = 0
        self.nsplit_right = 0
        self.mAddPercent = add_percent
        
    def __add__(self, other):
        self.nleft += other.nleft
        self.nright += other.nright
        self.noverlap += other.noverlap
        self.nunique_left += other.nunique_left
        self.nunique_right += other.nunique_right
        self.nidentical += other.nidentical
        self.nhalf += other.nhalf
        self.nsplit_left += other.nsplit_left
        self.nsplit_right += other.nsplit_right
        return self
    
    def getHeader( self ):
        h = "total_left\ttotal_right\tnoverlap\tnidentical\tnhalf\tunique_left\tunique_right\tsplit_left\tsplit_right"        
        if self.mAddPercent:
            h += "\t" + self.getHeaderPercent()
        return h

    def __str__( self ):
        h = "\t".join(map(str, \
                             (self.nleft, self.nright,
                              self.noverlap, self.nidentical, self.nhalf,
                              self.nunique_left, self.nunique_right,
                              self.nsplit_left, self.nsplit_right)))
        if self.mAddPercent:
            h += "\t" + self.asPercent()
            
        return h 

    def getHeaderPercent( self ):
        return "\t".join( map( lambda x: "pl%s\tpr%s" % (x,x), ("overlap", "identical", "half", "unique", "split")))
    
    def asPercent( self ):
        return "\t".join( map( lambda x: self.mPercentFormat % (100.0 * x),
                               ( float(self.noverlap) / self.nleft,
                                 float(self.noverlap) / self.nright,    
                                 float(self.nidentical) / self.nleft,
                                 float(self.nidentical) / self.nright,    
                                 float(self.nhalf) / self.nleft,
                                 float(self.nhalf) / self.nright,    
                                 float(self.nunique_left) / self.nleft,
                                 float(self.nunique_right) / self.nright,    
                                 float(self.nsplit_left) / self.nleft,
                                 float(self.nsplit_right) / self.nright )))

def getFile( options, section ):

    if options.output_filename_pattern:
        outfile = IOTools.openFile(options.output_filename_pattern % section, "w")
        E.info( "output for section '%s' goes to file %s" % (section, options.output_filename_pattern % section) )
    else:
        outfile = options.stdout
        outfile.write( "## section: %s\n" % section)
    return outfile

def writeDiff( outfile, symbol, genes ):

    for gene in genes:
        for exon in gene:
            outfile.write("%s\t%s\n" % (symbol, str(exon)))

def main( argv = None ):

    if not argv: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: diff_gtf.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-e", "--write-equivalent", dest="write_equivalent", 
                      help="write equivalent entries [default=%default].", action="store_true"  )

    parser.add_option("-f", "--write-full", dest="write_full",
                      help="write full gff entries [default=%default].", action="store_true"  )

    parser.add_option("-o", "--format=", dest="format",
                      help="output format [flat|multi-line] [default=%default]" )

    parser.add_option("-p", "--add-percent", dest="add_percent", action="store_true",
                      help="add percentage columns [default=%default]." )

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand", action="store_true",
                      help="ignore strand information [default=%default]." )

    parser.set_defaults(
        write_equivalent = False,
        write_full = False,
        format="flat",
        add_percent=False,
        ignore_strand = False,
        as_gtf = False,
        )

    (options, args) = E.Start( parser, argv, add_output_options = True )

    if len(args) != 2:
        raise ValueError( "two arguments required" )

    input_filename1, input_filename2 = args

    ## duplicated features cause a problem. Make sure
    ## features are non-overlapping by running
    ## gff_combine.py on GFF files first.

    E.info( "reading data started" )

    idx, genes2 = {}, set()
    for e in GTF.readFromFile( IOTools.openFile( input_filename2, "r" ) ):
        genes2.add( e.gene_id )
        if e.contig not in idx: idx[e.contig] = bx.intervals.intersection.Intersecter()
        idx[e.contig].add_interval( bx.intervals.Interval(e.start,e.end,value=e) )

    overlaps_genes = []

    E.info( "reading data finished: %i contigs" % len(idx))

    # outfile_diff and outfile_overlap not implemented
    # outfile_diff = getFile( options, "diff" )
    # outfile_overlap = getFile( options, "overlap" )
    overlapping_genes = set()

    genes1 = set()

    # iterate over exons
    with IOTools.openFile( input_filename1, "r" ) as infile:
        for this in GTF.iterator( infile ):

            genes1.add( this.gene_id )

            try:
                intervals = idx[this.contig].find( this.start, this.end )
            except KeyError:
                continue

            others = [ x.value for x in intervals ]
            for other in others:
                overlapping_genes.add( (this.gene_id, other.gene_id) )

            # check for identical/half-identical matches
            output = None
            for other in others:
                if this.start == other.start and this.end == other.end:
                    output, symbol = other, "="
                    break
            else:
                for other in others:
                    if this.start == other.start or this.end == other.end:
                        output, symbol = other, "|"
                        break
                else:
                    symbol = "~"

    # if outfile_diff != options.stdout: outfile_diff.close()
    # if outfile_overlap != options.stdout: outfile_overlap.close()

    outfile = None
    ##################################################################
    ##################################################################
    ##################################################################
    ## print gene based information
    ##################################################################    
    if overlapping_genes:
        outfile = getFile( options, "genes_ovl" )
        outfile.write ("gene_id1\tgene_id2\n" )
        for a, b in overlapping_genes: outfile.write( "%s\t%s\n" % (a,b) )
        if outfile != options.stdout: outfile.close()

        outfile_total = getFile(options, "genes_total" )
        outfile_total.write("set\tngenes\tnoverlapping\tpoverlapping\tnunique\tpunique\n" )

        outfile = getFile( options, "genes_uniq1" )
        b = set( [ x[0] for x in overlapping_genes ] )
        d = genes1.difference(b)
        outfile.write ("gene_id1\n" )
        outfile.write( "\n".join( d ) + "\n" )
        if outfile != options.stdout: outfile.close()
        outfile_total.write( "%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
                os.path.basename(input_filename1), len(genes1), len(b), 100.0 * len(b) / len(a),
                len(d), 100.0 * len(d) / len(genes1) ) )

        outfile = getFile( options, "genes_uniq2" )
        b = set( [ x[1] for x in overlapping_genes ] )
        d = genes2.difference(b)
        outfile.write ("gene_id2\n" )
        outfile.write( "\n".join( d ) + "\n" )
        if outfile != options.stdout: outfile.close()

        outfile_total.write( "%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
                os.path.basename(input_filename2), len(genes2), len(b), 100.0 * len(b) / len(a),
                len(d), 100.0 * len(d) / len(genes2) ) )
        if outfile_total != options.stdout: outfile_total.close()
        
    E.Stop()    

if __name__ == "__main__":
    sys.exit( main() )
