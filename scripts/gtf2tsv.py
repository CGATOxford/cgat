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
gtf2tsv.py - convert gtf file to a tab-separated table
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Genesets

Purpose
-------

convert gtf formatted file to tab-separated table with column headers for 
table import.

Note that coordinates are converted to 0-based open/closed notation (all on
the forward strand).

If -a/--attributes is set, attributes are converted into separate columns.

Usage
-----

Example::

   python gtf2tsv.py < in.gtf > out.tsv

Type::

   python gtf2tsv.py --help

for command line help.

Command line options
---------------------

'''
import os
import sys
import string
import re
import optparse
import CGAT.GTF as GTF
import CGAT.Experiment as E

def main():
    '''
    main function
    '''
    parser = E.OptionParser( version = "%prog version: $Id: gtf2tsv.py 2887 2010-04-07 08:48:04Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-o", "--only-attributes", dest="only_attributes", action="store_true",
                       help="output attributes as separate columns [default=%default]." )
    parser.add_option( "-f", "--full", dest="full", action="store_true",
                       help="output attributes as separate columns [default=%default]." )
    parser.add_option( "-i", "--invert", dest="invert", action = "store_true",
                       help="convert tab-separated table back to gtf [default=%default]." ) 
    parser.add_option( "-m", "--map", dest="map", type="choice",
                       choices=("transcript2gene", "peptide2gene", "peptide2transcript"),
                       help="output a map mapping transcripts to genes [default=%default]." ) 



    parser.set_defaults(
        only_attributes = False,
        full = False,
        invert = False,
        map = None,
        )


    (options,args) = E.Start( parser )

    if options.full:

        # output full table with column for each attribute
        attributes = set()
        data = []
        for gtf in GTF.iterator( options.stdin ):
            data.append(gtf)
            attributes = attributes.union( set(gtf.keys()) )

        # remove gene_id and transcript_id, as they are used
        # explicitely later
        attributes.difference_update( ["gene_id", "transcript_id"] ) 
            
        attributes = sorted(list(attributes))

        if options.only_attributes:
            header = [ "gene_id", "transcript_id"] + attributes
        else:
            header = [ "contig","source","feature","start","end","score","strand","frame","gene_id",
                       "transcript_id", ] + attributes

        options.stdout.write("\t".join( header ) + "\n" )

        if options.only_attributes:
            for gtf in data:
                options.stdout.write( "\t".join( map(str, (gtf.gene_id,
                                                           gtf.transcript_id,
                                                           ))))
                for a in attributes:
                    if a in ("gene_id", "transcript_id"): continue
                    try:
                        val = getattr( gtf, a )
                    except AttributeError:
                        val = ""
                    options.stdout.write( "\t%s" % val )
                options.stdout.write("\n")
        else:
            for gtf in data:
                options.stdout.write( "\t".join( map(str, (gtf.contig,
                                                           gtf.source,
                                                           gtf.feature,
                                                           gtf.start,
                                                           gtf.end,
                                                           gtf.score,
                                                           gtf.strand,
                                                           gtf.frame,
                                                           gtf.gene_id,
                                                           gtf.transcript_id,
                                                           ))))
                for a in attributes:
                    try:
                        val = getattr( gtf, a )
                    except AttributeError:
                        val = ""
                    options.stdout.write( "\t%s" % val )
                options.stdout.write("\n")
                
    elif options.invert:

        gtf = GTF.Entry()
        header = None
        for line in options.stdin:
            if line.startswith("#"): continue
            data = line[:-1].split("\t")
            if not header:
                header = data
                map_header2column = dict( [(y,x) for x,y in enumerate( header )] )
                continue
            
            # fill gtf entry with data
            try:
                gtf.contig = data[map_header2column["contig"]]
                gtf.source = data[map_header2column["source"]]
                gtf.feature = data[map_header2column["feature"]]
                # subtract -1 to start for 0-based coordinates
                gtf.start = int(data[map_header2column["start"]])
                gtf.end = int(data[map_header2column["end"]])
                gtf.score = data[map_header2column["score"]]
                gtf.strand = data[map_header2column["strand"]]
                gtf.frame = data[map_header2column["frame"]]
                gtf.gene_id = data[map_header2column["gene_id"]]
                gtf.transcript_id = data[map_header2column["transcript_id"]]
                gtf.parseInfo( data[map_header2column["attributes"]], line )
            except KeyError, msg:
                raise KeyError( "incomplete entry %s: %s: %s" % (str(data), str(map_header2column), msg))
            # output gtf entry in gtf format
            options.stdout.write( "%s\n" % str(gtf) )

    elif options.map:
        
        if options.map == "transcript2gene":
            fr = lambda x: x.transcript_id
            to = lambda x: x.gene_id
            options.stdout.write("transcript_id\tgene_id\n" )            
        elif options.map == "peptide2gene":
            fr = lambda x: x.protein_id
            to = lambda x: x.gene_id
            options.stdout.write("peptide_id\tgene_id\n" )            
        elif options.map == "peptide2transcript":
            fr = lambda x: x.protein_id
            to = lambda x: x.transcript_id
            options.stdout.write("peptide_id\ttranscript_id\n" )            
            
        map_fr2to = {}
        for gtf in GTF.iterator( options.stdin ):
            try:
                map_fr2to[fr(gtf)] = to(gtf)
            except AttributeError:
                pass

        for x,y in sorted(map_fr2to.iteritems()):
            options.stdout.write("%s\t%s\n" % (x,y) )
    else:
        header = ( "contig","source","feature","start","end","score","strand","frame","gene_id","transcript_id","attributes"  )
        options.stdout.write("\t".join( header ) + "\n" )

        for gtf in GTF.iterator( options.stdin ):
            
            attributes = []
            for a in gtf.keys():
                if a in ("gene_id", "transcript_id"): continue
                attributes.append( '%s %s' % (a, GTF.quote(gtf[a])) )

            attributes = "; ".join( attributes )
            
            options.stdout.write( "\t".join( map(str, (gtf.contig,
                                                       gtf.source,
                                                       gtf.feature,
                                                       gtf.start,
                                                       gtf.end,
                                                       GTF.toDot( gtf.score ),
                                                       gtf.strand,
                                                       gtf.frame,
                                                       gtf.gene_id,
                                                       gtf.transcript_id,
                                                       attributes,
                                                       ) ) ) + "\n" )
    E.Stop()

if __name__ == '__main__':
    sys.exit(main())


