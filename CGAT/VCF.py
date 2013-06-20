################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
VCF.py - Tools for working with VCF files
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The Variant Call Format (:term:`vcf`) is described
here:

http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0


Code
----

'''
import sys, gzip


class VCFEntry:
    def __init__(self, data, samples ):

        assert len(data) == len(samples) + 9
        self.contig, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format = \
            data[:9]
        
        self.genotypes = dict( zip(samples, data[9:]))
        self.order = samples

    def __str__(self):
        return "\t".join( map(str, (\
                    self.contig, self.pos, self.id, self.ref, self.alt, self.qual, 
                    self.filter, self.info, self.format,
                    "\t".join( [self.genotypes[x] for x in self.order ] ) ) ) )


class VCFFile:

    def __init__(self, infile):

        self.infile = infile
        self.format = {}
        self.info = {}
        self.fileformat = None

        while 1:
            line = self.infile.readline()
            
            if line.startswith("##"): 
                self.addMetaFromLine( line[2:-1] )
                continue
            elif line.startswith("#CHROM"):
                self.samples = line[:-1].split("\t")[9:]
                continue
            elif line.startswith("#"):
                continue
            break
        self.line = line

    def writeHeader( self, outfile, order = None ):
        outfile.write("##fileformat=%s\n" % self.fileformat)
        for key,values in self.format.iteritems():
            outfile.write("##FORMAT=%s,%s\n" % (key, ",".join(values)))
        for key,values in self.info.iteritems():
            outfile.write("##INFO=%s,%s\n" % (key, ",".join(values)))
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" )
        if order:
            assert len(order) == len(self.samples), \
                "number of samples do not match: %i != %i" % (len(order), len(self.samples))
            outfile.write("\t".join(order))
        else:
            outfile.write("\t".join(self.samples))
        outfile.write("\n")

    def __iter__( self ):
        return self

    def addMetaFromLine( self, line ):

        key, value = line.split("=",1)
        if key == "INFO": 
            data = value.split(",")
            self.info[data[0]] = data[1:]
        elif key == "FORMAT":
            data = value.split(",")
            self.format[data[0]] = data[1:]
        elif key=="fileformat":
            self.fileformat = value

    def next( self ):
        
        data = self.line[:-1].split("\t")
        self.line = self.infile.readline()
        if not self.line: raise StopIteration
        return VCFEntry(data, self.samples )
            
if __name__ == "__main__":
    
    inf = VCFFile(sys.stdin)
    
    for x in inf:
        print str(x)
