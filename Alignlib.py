"""Alignlib convenience functions.
"""
import sys, re, os, tempfile
#import alignlib

def writePairAlignment( row, col, map_row2col ):
    """wrapper around writePairAlignment."""
    
    tf = tempfile.TemporaryFile( "w+b" )
    alignlib.writePairAlignment( tf, row, col, map_row2col)
    tf.seek(0)
    r = tf.readlines()
    tf.close()
    return "".join(r)


def writeAlignataCompressed( map_row2col ):
    """wrapper around old syntax for writeAlignataCompressed."""

    out = tempfile.TemporaryFile( "w+b" )
    alignlib.writeAlignataCompressed( out, map_row2col )
    out.seek(0)
    r = out.readline()
    out.close()
    return r.split("\t")
