'''Tasks associated with genomic intervals.'''

import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import sqlite3
import cStringIO
import Experiment as E
import Pipeline as P
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed
# import Stats
import pysam
import numpy
import gzip
import fileinput

############################################################
############################################################
############################################################
## Pipeline configuration
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

############################################################
############################################################
############################################################
def getBedLocations(filename):
    '''return a list of regions as (contig,start,end) tuples
    from a bed file'''
    fh = open(filename,"r")
    region_list = []
    n_regions = 0

    for line in fh:
        if line.strip() != "" and line[0] !="#":
            fields = line.split("\t")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            region_list.append((contig,start,end))
            n_regions +=1

    fh.close()

    #E.info("Read in %i regions from %s" % ( n_regions, filename) )
    return (region_list)

############################################################
############################################################
############################################################
def mergeBedFiles( infiles, outfile ):
    '''generic method for merging bed files. '''

    if len(infiles) < 2:
        raise ValueError( "expected at least two files to merge into %s" % outfile )

    infile = " ".join( infiles )
    statement = '''cat %(infile)s |
                   mergeBed -i stdin |
                   cut -f 1-3 |
                   awk '{printf("%%s\\t%%i\\n",$0, ++a); }'
                   > %(outfile)s ''' 
    P.run()

############################################################
############################################################
############################################################
def mergeIntervalsWithScores( infile, outfile, dist, method ):
    '''merge adjacent intervals (within dist) and integrate scores (using method) from merged peaks.
       Assume bed file sorted by position and score in column 5.
       Methods: mean, max, length weighted mean'''

    intervals = open( infile, "r")
    merged = open( outfile, "w")

    topline = intervals.readline()
    last_contig, last_start, last_end, last_id, last_score = topline[:-1].split("\t")[:5]
    last_start = int(last_start)
    last_end = int(last_end)
    last_score = int(last_score)
    for line in intervals:
        data = line[:-1].split("\t")
        contig, start, end, interval_id, score = data[:5]
        start = int(start)
        end = int(end)
        score = int(score)
        if (contig == last_contig) and ((last_end+dist) >= start):
            if method == "mean":
                newscore = (score + last_score) / 2
            elif method == "length_weighted_mean":
                length1 = end - start
                length2 = last_end - last_start
                newscore = ((score*length1)+(last_score*length2))/(length1+length2)
            elif method == "max":
                newscore = max(last_score, score)
            last_end=end
            last_score=newscore
        else:
            merged.write("%(last_contig)s\t%(last_start)i\t%(last_end)i\t%(last_id)s\t%(last_score)s\n" % locals())
            data = line[:-1].split("\t")
            last_contig, last_start, last_end, last_id, last_score = data[:5]
            last_start = int(last_start)
            last_end = int(last_end)
            last_score = int(last_score)
    intervals.close()
    merged.close()

############################################################
############################################################
############################################################
def intersectBedFiles( infiles, outfile ):
    '''generic method for merging bed files.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered.
    '''

    if len(infiles) == 1:
        shutil.copyfile( infiles[0], outfile )

    elif len(infiles) == 2:
        
        statement = '''
        intersectBed -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % (infiles[0], infiles[1])

        P.run()
        
    else:

        tmpfile = P.getTempFilename(".")

        # need to merge incrementally
        fn = infiles[0]
        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()
        
        for fn in infiles[1:]:
            statement = '''mergeBed -i %(fn)s | intersectBed -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %(outfile)s '''
        P.run()

        os.unlink( tmpfile )

############################################################
############################################################
############################################################
def subtractBedFiles( infile, subtractfile, outfile ):
    '''subtract intervals in *subtractfile* from *infile*
    and store in *outfile*.
    '''
    statement = '''
        intersectBed -v -a %(infile)s -b %(subtractfile)s |\
        cut -f 1,2,3,4,5 |\
        awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % locals()

    P.run()

############################################################
############################################################
############################################################
def makeIntervalCorrelation( infiles, outfile, field, reference ):
    '''compute correlation of interval properties between sets  '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip( infile, ".bed" )
        tablename = "%s_intervals" % P.quote( track )
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals()
        cc.execute( statement )
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add( contig, start, end, peakval )        
        idx.append( ix )
        tracks.append( track )
    outs = open( outfile, "w" )
    outs.write( "contig\tstart\tend\tid\t" + "\t".join( tracks ) + "\n" )

    for bed in Bed.iterator( infile = open( reference, "r") ):
        
        row = []
        for ix in idx:
            try:
                intervals = list(ix.get( bed.contig, bed.start, bed.end ))
            except KeyError:
                row.append( "" )
                continue
        
            if len(intervals) == 0:
                peakval = ""
            else:
                peakval = str( (max( [ x[2] for x in intervals ] )) )
            row.append( peakval )

        outs.write( str(bed) + "\t" + "\t".join( row ) + "\n" )

    outs.close()

