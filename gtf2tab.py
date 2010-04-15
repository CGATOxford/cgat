################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: gtf2tab.py 2887 2010-04-07 08:48:04Z andreas $
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
import os, sys, string, re, optparse

USAGE="""python %s [OPTIONS] < in.gtf

convert gtf formatted file to tab-separated table with column headers for 
table import.

Note that coordinates are converted to 0-based open/closed notation (all on
the forward strand).

If -a/--attributes is set, attributes are converted into separate columns.

"""

import GFF, GTF
import Experiment as E

def main():
    '''
    main function
    '''
    parser = optparse.OptionParser( version = "%prog version: $Id: gtf2tab.py 2887 2010-04-07 08:48:04Z andreas $", usage = USAGE)

    parser.add_option( "-o", "--only-attributes", dest="only_attributes", action="store_true",
                       help="output attributes as separate columns [default=%default]." )
    parser.add_option( "-f", "--full", dest="full", action="store_true",
                       help="output attributes as separate columns [default=%default]." )
    parser.add_option( "-i", "--invert", dest="invert", action = "store_true",
                       help="convert tab-separated table back to gtf [default=%default]." ) 

    parser.set_defaults(
        only_attributes = False,
        full = False,
        invert = False,
        )


    (options,args) = E.Start( parser )

    if options.full:
        # output full table with column for each attribute
        attributes = set()
        data = []
        for gtf in GTF.iterator( options.stdin ):
            data.append(gtf)
            attributes = attributes.union( set(gtf.keys()) )
            
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

    else:
        header = ( "contig","source","feature","start","end","score","strand","frame","gene_id","transcript_id","attributes"  )
        options.stdout.write("\t".join( header ) + "\n" )

        for gtf in GTF.iterator( options.stdin ):
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
                                                       str(gtf.attributes),
                                                       ) ) ) + "\n" )
            E.Stop()

if __name__ == '__main__':
    sys.exit(main())


