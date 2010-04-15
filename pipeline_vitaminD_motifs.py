import sys, re, os, tempfile, collections

import Experiment as E
import Pipeline as P
import sqlite3

PARAMS = P.getParameters()

def filterMotifsFromMEME( infile, outfile, selected ):
    '''select motifs from a MEME file and save into outfile
    '''

    outs = open(outfile, "w" )
    if len(selected) == 0:
        outs.close()
        return

    keep = True

    for line in open(infile, "r"):
        if line.startswith( "MOTIF" ):
            motif = re.match( "MOTIF\s+(\d+)", line ).groups()[0]
            if motif in selected:
                keep = True
            else:
                keep = False
        if keep: outs.write( line )
    outs.close()
