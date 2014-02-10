'''
transfac2transfac.py - filter transfac motif files
====================================================

:Author: Stephen Sansom
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Filter a transfac motif file.

Usage
-----

Example::

   python cgat_script_template.py 

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

# parse water alignments to summary stats in a tab delimited text format.

# this is what we are trying to parse

##AC
##
##ID V$.....
##
##//

import sys, os, re
import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-f", "--filter", dest="filter", default=False,
                      help="ID prefix to filter on, eg. V for vertebrates")

    parser.add_option("-p", "--pattern", dest="pattern", default=False,
                      help="ID pattern to filter (filter is case insensitive) eg. pax6. Multiple patterns should be specified as a comma separated list")

    (options,args) = E.Start( parser )

    if options.pattern != False:
        patterns = [ x.strip() for x in options.pattern.split(",") ]
        E.info("Supplied patterns %s" % ", ".join(patterns))
    else:
        patterns = False

    filtered_motifs=[]
    n=0

    inmotif, tid, filter_emit, pattern_emit = False, False, False, False

    for line in options.stdin:

        # pick up motif start and ends.
        if line.startswith("AC") and inmotif == False:
            #print "in align"
            inmotif = True
            motif = line
            continue
        elif line.startswith("ID") and inmotif == True:
            #print line
            tid = line.split("  ")[1]
            motif += line
            continue

        elif line.startswith("//") and inmotif == True:

            motif += line

            if tid == False:
                raise ValueError("matrix ID not determined")


            if options.filter != False:
                if tid.startswith(options.filter):
                    filter_emit = True
            else:
                filter_emit = True

            if patterns != False:
                for pat in patterns:
                    match = re.search(pat,tid,re.IGNORECASE)
                    if match != None:
                        pattern_emit = True
                        break
            else:
                pattern_emit = True


            if filter_emit==True and pattern_emit==True:
                filtered_motifs.append(motif)
                n+=1

            inmotif, tid, filter_emit, pattern_emit = False, False, False, False
            continue

        elif inmotif == True:
            motif += line

        elif inmotif == False:
            continue

        else:
            raise ValueError("unknown parsing state")
        
    options.stdout.write("".join(filtered_motifs))

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

