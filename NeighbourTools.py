####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: NeighbourTools.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


import math
import Tools
import numpy

import alignlib

from Table_nrdb import Table_nrdb
from TablePairsdbNeighbours import TablePairsdbNeighbours

def GetAllRanges( tbl_pairsdb_90x90, nids, lengths ):
    """retrieve all overlaps between a list of neighbours.
    """

    nsequences = len(nids)
    all_ranges = []
    
    for x in range(0, nsequences):
        all_ranges.append([])
        nid1 = nids[x]
        
        for y in range(0, nsequences):
            nid2 = nids[y]

            if nid1 == nid2:
                all_ranges[x].append( [(1, lengths[x]),] )
            else:
                xranges = Tools.CombineIntervallsLarge(list(tbl_pairsdb_90x90.GetBoundaries( nid1, nid2 )))
                if not xranges:
                    all_ranges[x].append( [] )
                else:
                    all_ranges[x].append( xranges )

    return all_ranges

def GetAllLengths( tbl_nrdb, nids ):
    """retrieve all lengths for a set of nids.
    """
    lengths = []
    for nid in nids:
        lengths.append( tbl_nrdb.GetLength( nid ))

    return lengths

##-----------------------------------------------------------------------------------------------------
def BuildNeighbourBLASTMatrix( dbhandle, query_nid, resolution = 1.0, table_name = None):
    """build neighbourhood-matrix based on BLAST alignments
    to query_nid.
    """

    tbl_nrdb = Table_nrdb(dbhandle)
    tbl_pairsdb_90x90 = TablePairsdbNeighbours( dbhandle )
    if table_name:
        tbl_pairsdb_90x90.SetName( table_name )

    nids = [query_nid,]
    nids += map(lambda x: x[3], tbl_pairsdb_90x90.GetNeighbours( query_nid, sort_order = 3, skip_query = 1 ))
    
    nsequences = len(nids)
    
    ##---------------------------------------------------------------------------
    ## get sequence lengths
    lengths = NeighbourTools.GetAllLengths( tbl_nrdb, nids )

    x = 0
    start_points = []
    for l in lengths:
        start_points.append(x/resolution)
        x+= l
    total_length = x
    
    ##---------------------------------------------------------------------------        
    # retrieve sequence ranges
    all_ranges = GetAllRanges( tbl_pairsdb_90x90, nids, lengths )
        
    matrix = numpy.zeros( (nsequences, total_length / resolution), numpy.UnsignedInt8)    

    current_row = 0
    for y in range(0, nsequences):

        for x in range(0,nsequences):
            for r_from,r_to in all_ranges[x][y]:
                r_from /= param_resolution
                r_to /= param_resolution                
                matrix[current_row,(start_points[x] + r_from-1): start_points[x] + r_to] = 1.0
                
        current_row += 1
    

    return start_points, matrix

##-----------------------------------------------------------------------------------------------------
def BuildBLASTMatrix( dbhandle,
                      query_nid,
                      resolution = 1.0,
                      table_name = None,
                      combine_repeats = None,
                      max_evalue = None,
                      min_evalue = None,
                      residue_level = None,
                      parser = None,
                      add_self = None):
    """build matrix based on BLAST alignments to query_nid.

    matrix of size N*M
    N: number of neighbours
    M: length of query (scaled with resolution)

    alignments are truncated.

    the query is included in the matrix.
    
    if combine_repeats is set, multiple alignments between the query and a sbjct will
    be entered into the same row.

    if residue_level is set, entries are added on the residue level. The resolution parameter
    is ignored.
    """

    if residue_level:
        query_length = Table_nrdb(dbhandle).GetLength( query_nid )
    else:
        query_length = int( math.floor( float(Table_nrdb(dbhandle).GetLength( query_nid )) / float(resolution)))

    tbl_pairsdb_90x90 = TablePairsdbNeighbours( dbhandle )
    if table_name:
        tbl_pairsdb_90x90.SetName( table_name )

    neighbours = tbl_pairsdb_90x90.GetNeighbours( query_nid,
                                                  sort_order = 3,
                                                  skip_query = add_self,
                                                  min_evalue = min_evalue,
                                                  max_evalue = max_evalue)
    nindex = {}
    
    nneighbours = 0
    if combine_repeats:
        for neighbour in neighbours:
            (query_from, query_to, query_ali,
             sbjct_nid, sbjct_from, sbjct_to, sbjct_ali, score, pide, evalue) = neighbour
            if not nindex.has_key(sbjct_nid):
                nindex[sbjct_nid] = nneighbours
                nneighbours += 1
    else:
        nneighbours = len(neighbours)

    if add_self:
        nneighbours += 1

    matrix = numpy.zeros( (nneighbours, query_length), numpy.int)    

    if add_self:
        matrix[0, 0:query_length] = 1
        row = 1
    else:
        row = 0
        
    for neighbour in neighbours:

        (query_from, query_to, query_ali,
         sbjct_nid, sbjct_from, sbjct_to, sbjct_ali, score, pide, evalue) = neighbour

        if combine_repeats:
            use_row = nindex[sbjct_nid]
        else:
            use_row = row
            row += 1

        if residue_level:
            map_sbjct2query = alignlib.makeAlignataVector()
            alignlib.fillAlignataCompressed( map_sbjct2query, sbjct_from, sbjct_ali, query_from, query_ali )
            if parser:
                parser( map_sbjct2query )
                
            for x in range(sbjct_from, sbjct_to + 1):
                y = map_sbjct2query.mapRowToCol(x)
                if y:
                    try:
                        matrix[use_row, y-1] = 1
                    except IndexError:
                        print "IndexError in ", query_nid, sbjct_nid, x, y-1, query_length
        else:
            yfrom = int(math.floor(query_from/resolution))
            yto   = int(math.floor(query_to/resolution)) 
            matrix[use_row, yfrom:yto] = 1
            
    return matrix

##-----------------------------------------------------------------------------------------------------
def GetNeighbourMatrixNids( dbhandle,
                            nids,
                            table_name = None,
                            max_evalue = None,
                            use_value = None,
                            method = "evalue"):
    """build matrix between proteins based based on BLAST hits.

    matrix of size N*N
    N: number of sequences

    if combine_repeats is set to 1, multiple alignments between the query and a sbjct will
    be entered into the same row.
    """

    nnids = len(nids)
    rows = {}
    for nid in nids:
        rows[nid] = len(rows)
        
    tbl_pairsdb_90x90 = TablePairsdbNeighbours( dbhandle )
    if table_name:
        tbl_pairsdb_90x90.SetName( table_name )

    matrix = numpy.zeros( (nnids, nnids), numpy.float)
    
    map_row2nid = []
    
    for nid in nids:

        neighbours = tbl_pairsdb_90x90.GetNeighboursStrength( nid, 
                                                              max_evalue = max_evalue,
                                                              select_statement = method)
        
        for neighbour in neighbours:

            (sbjct_nid, evalue) = neighbour
            
            if not rows.has_key(sbjct_nid):
                continue
            
            if use_value:
                matrix[rows[nid],rows[sbjct_nid]] = use_value
            else:
                matrix[rows[nid],rows[sbjct_nid]] = evalue


    
    return matrix, rows




