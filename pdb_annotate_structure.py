################################################################################
#   Gene prediction pipeline 
#
#   $Id: pdb_annotate_structure.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
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


USAGE="""
python pdb_annotate_structure.py [OPTIONS] pdb_id

annotate a structure given by pdb_id
"""


import sys, string, re, optparse

import alignlib

import Experiment
import IOTools
import PdbTools 

DEBUG = 1

#--------------------------------------------------------------------------------
def MapRight( mapping, row_residue ):
    """return mapping of row_residue Go right, if not found."""
    col_residue = mapping.mapColToRow( row_residue )
    
    max_residue = mapping.getColTo()

    while col_residue == 0:
        row_residue = row_residue + 1
        if row_residue > max_residue:
            return (0)
            
        col_residue = mapping.mapColToRow( row_residue )

    return (col_residue)

#--------------------------------------------------------------------------------
def MapLeft( mapping, row_residue ):
    """return mapping of row_residue Go right, if not found."""
    col_residue = mapping.mapColToRow( row_residue )
    
    min_residue = mapping.getColFrom()
    
    while col_residue == 0:
        row_residue = row_residue - 1
        if row_residue < min_residue:
            return (0)
            
        col_residue = mapping.mapColToRow( row_residue )

    return (col_residue)


if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: pdb_annotate_structure.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = USAGE )

    parser.add_option( "-p", "--filename-pdb", dest="filename_pdb", type="string",
                       help="filename with pdb structure."  )

    parser.add_option( "-f", "--filename-fasta", dest="filename_fasta", type="string",
                       help="filename with sequence on which annotation is based. If not given, the pdb sequence is used."  )

    parser.set_defaults(
        filename_pdb = None,
        filename_fasta = None,
        )
    
    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if len(args) != 1:
        print USAGE, "please supply the pdb identifier"
        sys.exit(2)

    if string.find(args[0], "-") != -1:
        (options.pdb_id,options.pdb_chain) = string.split(args[0], "-")
    else:
        options.pdb_id = args[0]
        options.pdb_chain = ""

    message = ""
    
    ## pdb_filename = "/homes/heger/pdb_temp/temp.pdb"
    ## print PdbTools.GetPdbFile( pdb_id, pdb_filename )
    ## pdb_lines = PdbTools.GetPdbFileLine( pdb_id )

    #######################################################################
    ## retrieve structure
    if options.filename_pdb:
        infile = open(options.filename_pdb, "r")
        pdb_lines = infile.readlines()
        infile.close()
    else:
        pdb_lines = os.popen( param_retrieval_command % string.lower(param_pdb_id) ).readlines()
    
    viewer = PdbTools.RasmolViewInline( pdb_lines, sys.stdout )
    viewer.Command( "echo %s" % message)

    if options.filename_fasta:
        infile = open( options.filename_fasta, "r" )
        description, reference_sequence = IOTools.readSequence( infile )
        infile.close()
    else:
        reference_sequence = None
    
    if DEBUG:
        viewer.Command("echo cmdline: %s" % (string.join(sys.argv, " ")))
                      
    if not pdb_lines:
        viewer.Command("echo error: structure not found in local database")
        viewer.WriteScript()
        sys.exit()


    if reference_sequence:
        map_pdb2seq, rmap_pdb2seq, rmap_seq2pdb, lstructure, first_residue, last_residue, sequence = PdbTools.buildMapPdb2Sequence( reference_sequence,
                                                                                                                                    options.filename_pdb,
                                                                                                                                    options,
                                                                                                                                    options.pdb_chain )
        ## note: switchRowCol does not work for AlignataVector, thus the inefficient
        ## mapping of col to row.
    else:
        map_pdb2seq, map_pdb2seq = None, None

    first = True
    ninput, noutput, nmissed = 0, 0, 0
    
    for line in sys.stdin:
        if line[0] == "#" : continue
        if first:
            first = False
            continue

        residue_number, aminoacid, shape, colour = line[:-1].split("\t")

        if colour == "": colour = None
        if shape == "": shape = None
        
        ninput += 1

        # use pdb2seq for mapping, as the script renumbers the residues.
        if "-" in residue_number:
            first_res, last_res = residue_number.split( "-" )
            pdb_from = MapRight( map_pdb2seq, int(first_res) )
            pdb_to = MapLeft( map_pdb2seq, int(last_res) )
        else:
            pdb_from = map_pdb2seq.mapColToRow( int(residue_number) )
            pdb_to = pdb_from

        if pdb_from == 0 or pdb_to == 0:
            nmissed += 1
            if options.loglevel >= 1:
                options.stdlog.write( "# residue not found: %s:%s\t%i-%i\n" % (residue_number, aminoacid, pdb_from, pdb_to ))
            continue
        
        if options.loglevel >= 3:
            options.stdlog.write( "# mapped: %s%s to %s%i-%s%i\n" % (aminoacid, residue_number, sequence[pdb_from-1], pdb_from, sequence[pdb_to-1], pdb_to) )

        viewer.highlightResidues( int(pdb_from), int(pdb_to), options.pdb_chain,
                                  colour = colour,
                                  shape = shape )

        noutput += 1
            
    if options.pdb_chain:
        viewer.Command( "restrict :%s" % options.pdb_chain)
        viewer.Command( "centre" )

    viewer.WriteScript()

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nmissed=%i\n" % (ninput, noutput, nmissed))
        
    Experiment.Stop()
    






