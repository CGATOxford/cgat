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
gff2exons.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gff2exons.py --help

Type::

   python gff2exons.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse

USAGE="""python %s [OPTIONS] < in.gff > reference.exons

Convert GFF list to exons. 
Can also create a file peptides2genes.

Version: $Id: gff2exons.py 2447 2009-01-27 17:12:48Z andreas $
""" % sys.argv[0]

import CGAT.Experiment as Experiment
import CGAT.GFF as GFF
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta

def processEntries( name, entries, options, fasta, contigs ):

    ## reorder, if negative strand
    # if Genomics.IsNegativeStrand( entries[0].strand ):
    # entries.reverse()

    is_negative = Genomics.IsNegativeStrand(entries[0].strand)

    contig = entries[0].contig

    lcontig = contigs[contig]

    # sort in-order in transcript
    entries.sort( key=lambda x: x.start )
    if is_negative: entries.reverse()

    for gff in entries:
        if gff.end > lcontig or gff.start >= lcontig:
            if options.loglevel >= 1:
                options.stdlog.write( "# coordinates for %s on %s out of bounds (%i:%i > %i)\n" % \
                                          (str(gff.mAttributes), contig, gff.start, gff.end, lcontig) )
            return False

        
    if options.convert_to_cds:
        cds_start = 0
        t = 0
        for gff in entries:
            t += gff.end - gff.start
            
        options.stdout.write( "\t".join( map(str, (
            name,
            entries[0].contig,
            "+",
            1,
            0,
            cds_start, t,
            cds_start, t,
            ))) + "\n" )
    else:
        n = 0
        cds_start = 0
        cds_end = 0

        if options.reset_coordinates:
            if is_negative:
                offset = -(entries[-1].start)
            else:
                offset = -(entries[0].start)
        else:
            offset = 0

        for gff in entries:
            n += 1
            cds_end += gff.end - gff.start

            if offset:
                gff.start += offset
                gff.end += offset

            options.stdout.write( "\t".join( map(str, (
                name,
                gff.contig,
                gff.strand,
                gff.frame,
                n,
                cds_start, cds_end,
                gff.start,
                gff.end ))) + "\n" )

            cds_start = cds_end

    return True

##------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gff2exons.py 2447 2009-01-27 17:12:48Z andreas $")

    parser.add_option("-p", "--filename-peptides2genes", dest="filename_peptides2genes", type="string",
                      help="filename in which to output peptides2genes information."  )
    parser.add_option("-r", "--remove-unknown-contigs", dest="remove_unknown_contigs", action="store_true",
                      help="remove contigs that are in the genomic database."  )
    parser.add_option("-s", "--reset-coordinates", dest="reset_coordinates", action="store_true",
                      help="all coordinates start 0. Use this for cds based assemblies."  )
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genomic data (indexed)." )
    parser.add_option( "--convert-to-cds", dest="convert_to_cds", action="store_true",
                      help="convert to cds coordinates. Introns are removed and all predictions are located onto positive strand."  )
    
    parser.set_defaults(
        filename_peptides2genes = None,
        separator = " ",
        translate_contigs = True,
        remove_unknown_contigs = False,
        reset_coordinates = False,
        convert_to_cds = False,
        genome_file = None )

    (options, args) = Experiment.Start( parser )

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file ) 
        contigs = fasta.getContigSizes()
    else:
        fasta = None
        contigs = None

    if (options.remove_unknown_contigs or options.translate_contigs) and not fasta:
        raise ValueError("please supply genomic sequence.")

    gffs = GFF.iterator( sys.stdin )

    map_peptides2genes = {}

    ninput, noutput, nskipped, ntranslated = 0, 0, 0, 0

    last_name = None
    for gff in gffs:

        if gff.feature != "CDS": continue

        ninput += 1

        if contigs:
            if options.translate_contigs:
                if gff.contig not in fasta.mIndex and gff.contig in contigs:
                    gff.contig = fasta.mSynonyms[gff.contig]
                    ntranslated += 1

            if gff.contig not in contigs:
                if options.remove_unknown_contigs:
                    if options.loglevel >= 1:
                        options.stdlog.write( "# contig %s for %s unknown\n" %\
                                          (gff.contig, str(gff.mAttributes) ) )
                    
                    nskipped += 1
                    continue
                else:
                    raise IndexError("unknown contig %s in %s" % (gff.contig, str(gff)))

        name = gff.mAttributes["protein_id"]

        if name != last_name:
            if last_name:
                if processEntries( last_name, entries, options, fasta, contigs ):
                    noutput += 1
                else:
                    nskipped += 1

            entries = []

            if name in map_peptides2genes:
                raise ValueError( "file is not sorted correctly: %s occurs twice" % (name) )

            map_peptides2genes[ name ] = gff.mAttributes['gene_id']
        entries.append( gff )
        last_name = name

    if processEntries( last_name, entries, options, fasta, contigs ):
        noutput += 1
    else:
        nskipped += 1
        
    if options.filename_peptides2genes:
        outfile = open(options.filename_peptides2genes, "w")
        for peptide, gene in map_peptides2genes.items():
            outfile.write( "%s\t%s\n" % (peptide,gene) )
        outfile.close()
        
    genes = set(map_peptides2genes.values())
    
    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i, ntranslated=%i, ngenes=%i, npeptides=%i\n" %\
                                 (ninput,
                                  noutput,
                                  nskipped,
                                  ntranslated,
                                  len(genes),
                                  len(map_peptides2genes) ) )
    Experiment.Stop()



