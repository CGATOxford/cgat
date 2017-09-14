import os
import re

from SphinxReport.Renderer import *
from SphinxReport.Tracker import *

# for trackers_derived_sets and trackers_master
if not os.path.exists("conf.py"):
    raise IOError( "could not find conf.py" )

execfile( "conf.py" )

def iterate_with_wrap( lines ):
    """iterate over lines in a Makefile

    The iterator merges multi-line lines.
    Comments before a variable are collected.
    """
    
    comments, buf = [], []

    keep = False

    for line in lines:
        if re.search( "Section parameters: start", line ): 
            keep = True
            comments = []
            continue
        if re.search( "Section parameters: end", line ): 
            break

        if line.startswith("#"):
            comments.append( re.sub( "^#+", "", line[:-1] ) )
            continue

        if not line.strip():
            comments = []
            continue

        if not keep: continue

        if line.endswith("\\\n"):
            buf.append( line[:-2].strip() )
            continue

        if buf:
            buf.append( line.strip() )
            yield " ".join( buf ), "\n".join(comments)
            buf = []            
            comments = []

        else:
            yield line, "\n".join(comments)
            comments = []

class MakefileParameters( Tracker ):
    """list parameters in the Makefile of a pipeline."""

    @returnLabeledData
    def __call__( self, track, slice = None ):

        rx = re.compile("^(\S[^=]*)\s*=\s*(.*)")
        infile = open( os.path.join( data_dir, track), "r" )
        data = []
        for line, comments in iterate_with_wrap( infile ):

            a = rx.search( line )
            if not a: continue
            
            k,v = a.groups() 
            
            if k.endswith("?"): k=k[:-1]
            if not v: v = "<empty>"
            if comments:
                data.append( (k, v + "\n\n" + comments ) )
            else:
                data.append( (k, v) )

        infile.close()
        return data
