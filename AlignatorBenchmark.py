####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: AlignatorBenchmark.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


#--------------------------------------------------------------------------------------
# compare two clusterings.
#
# 
#
#

import sys
import re
import string
import os
import time 

from Picasso import *

import alignlib

from Table_nrdb import Table_nrdb
from Table_benchmark_sources import Table_benchmark_sources
from Table_benchmark_alignments import Table_benchmark_alignments

class AlignatorBenchmark:

    def __init__ (self, dbhandle ):
        self.dbhandle = dbhandle

        self.mLogLevel = 2
        
        # minimum overlap for two regions of the two clusters to overlap
        self.mMinOverlap = 5

        # where to store temporary results
        self.mTempFilename = PATH_TEMP + "/alignator_benchmark.tmp"

        self.mTableNameSources = None
        self.mTableNameAlignments = None

        self.mTableSources      = None
        self.mTableAlignments   = None

        self.mAlignator         = None

    #-------------------------------------------------------------------------------------------------------
    def SetTableSources( self, table_name ):
        self.mTableNameSources = table_name
        self.mTableSources = Table_benchmark_sources( self.dbhandle, self.mTableNameSources )
        self.mTableSources.Create()
        
    #-------------------------------------------------------------------------------------------------------
    def SetAlignator( self, alignator ):
        self.mAlignator = alignator
    
    #-------------------------------------------------------------------------------------------------------
    def FillSources( self, statement ):
        self.mTableSources.Fill( statement )

    #-------------------------------------------------------------------------------------------------------
    def FilterSources( self, where_statement):
        """remove entries of families, which appear more or less than a certain number
        of times.
        """
        self.mTableSources.DeleteSelectedClasses( where_statement )
        
    #-------------------------------------------------------------------------------------------------------
    def FilterSourcesDegeneracy( self, min_degeneracy, max_degeneracy ):
        """remove entries of families, which appear more or less than a certain number
        of times.
        """
        degeneracies = self.mTableSources.GetClassDegeneracies()

        for class_name, degeneracy in degeneracies:
            if degeneracy > max_degeneracy:
                if self.mLogLevel >= 2:
                    print "--> maximum exceeded: deleting %i from %s" % (degeneracy - max_degeneracy, class_name)
                self.mTableSources.DeleteSomeFromClass( class_name, degeneracy - max_degeneracy )
            elif degeneracy < min_degeneracy:
                if self.mLogLevel >= 2:
                    print "--> minimum exceeded: deleting %i from %s" % (degeneracy, class_name)
                self.mTableSources.DeleteSomeFromClass( class_name, degeneracy )
                
                
    #-------------------------------------------------------------------------------------------------------
    def SetTableAlignments( self, table_name ):
        self.mTableNameAlignments = table_name
        self.mTableAlignments = Table_benchmark_alignments( self.dbhandle, self.mTableNameAlignments )
        self.mTableAlignments.Create()

    #-------------------------------------------------------------------------------------------------------
    def FillAlignments( self ):
        """the main all-vs-all alignments engine.
        """
        sources = self.mTableSources.GetSources()

        nsources = len(sources)
        if self.mLogLevel >= 1:
            print "--> calculating %i alignments for %i sequences" % ((nsources * (nsources -1)) / 2, nsources)

        # calculate the alignanda objects
        if self.mLogLevel >= 1:
            print "--> retrieving alignanda objects"
            sys.stdout.flush()
            
        alignanda = self.CreateAlignandumObjects( sources )

        map_query2sbjct = alignlib.makeAlignataVector()
        outfile = open (self.mTempFilename, "w" )
        
        
        # do all vs all alignments
        for query in range(0, nsources - 1):
            query_id, query_alignandum = alignanda[query]

            start_time = time.time()
            
            if self.mLogLevel >= 2:
                print "processing id %i at %s" % (query_id, time.asctime(time.localtime(start_time)))
                sys.stdout.flush()
            
            for sbjct in range( query + 1, nsources):
                sbjct_id, sbjct_alignandum = alignanda[sbjct]                

                self.mAlignator.Align( query_alignandum, sbjct_alignandum, map_query2sbjct )

                (query_ali, sbjct_ali) = alignlib.writeAlignataCompressed( map_query2sbjct )

                outfile.write( string.join( map( str, (
                    query_id,
                    map_query2sbjct.getRowFrom(),
                    map_query2sbjct.getRowTo(),
                    query_ali,
                    sbjct_id,
                    map_query2sbjct.getColFrom(),
                    map_query2sbjct.getColTo(),
                    sbjct_ali,
                    map_query2sbjct.getScore(),
                    map_query2sbjct.getNumGaps(),
                    map_query2sbjct.getLength(),
                    0)), "@") + "\n" )
                    
                    
            stop_time = time.time()
            
            if self.mLogLevel >= 2:
                print "--> alignments: %5i, time: %7.2fs" %\
                      ( nsources - query - 1,
                        stop_time - start_time)

        outfile.close()

        # load data
        self.mTableAlignments.Drop()
        self.mTableAlignments.Create()
        self.mTableAlignments.Load( self.mTempFilename)        

        
    #-------------------------------------------------------------------------------------------------------
    # methods to overload
    #-------------------------------------------------------------------------------------------------------
    def CreateAlignandumObjects( self, sources ):

        tbl_nrdb = Table_nrdb( self.dbhandle )

        alignanda = []
        
        for id, nid, nid_from, nid_to in sources:
            
            if self.mLogLevel >= 2:
                print id,
                sys.stdout.flush()
                
            sequence = tbl_nrdb.Get_Sequence_From_NID( nid )
            alignandum = alignlib.makeSequence( sequence[nid_from-1:nid_to] )
            alignanda.append( (id, alignandum) )

        if self.mLogLevel >= 2:
            print
            
        return alignanda
    
    #-------------------------------------------------------------------------------------------------------        
        
if __name__ == '__main__':

    dbhandle = Picasso()
    if not dbhandle.Connect():
	print "Connection failed"
	sys.exit(1)

    x = AlignatorBenchmark( dbhandle )

    x.SetTableSources( "picasso_benchmark.scop_sources" )
    x.SetTableAlignments( "picasso_benchmark.scop_alignments" )

    # if a sequence has repeated domains of the same type, only one gets chosen
    statement = "SELECT DISTINCTROW rep_nid, rep_from, rep_to, SUBSTRING(scop_class, 1, 9) " +\
                " FROM nrdb40_scop GROUP BY rep_nid, SUBSTRING(scop_class, 1, 9) "
    
    x.FillSources( statement )

    x.FilterSourcesDegeneracy( 2, 3 )
    x.FilterSources( "SUBSTRING(class,1,3) > 4")

    alignator = alignlib.makeFullDP( -12.0, -2.0 )

    x.SetAlignator( alignator )
    x.FillAlignments()
    





