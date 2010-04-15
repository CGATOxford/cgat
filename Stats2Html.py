####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Stats2Html.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


# translate an ascii-table into an html-table

from Pairsdb import *

import sys
import os
import string
import HTMLgen
import re

from Table_nrdb import Table_nrdb
import Selects_1

from RSDB2HTML import RSDB2HTML

PICASSO_URL = "http://jura.ebi.ac.uk:8801/PairsDB"

def Build_report_sequence ( nid ):
    return "<A HREF='%s/report_sequence?nid=%s' target='_blank'>%s</A>" % (PICASSO_URL, nid, nid )

def Build_scop ( scop_id ):
    domain_id = scop_id + "_" * (6 - len(scop_id))
    return "<A HREF='http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?sid=d%s' target='_blank'>%s</A>" % (domain_id, scop_id )

def Build_join( nid1, nid2 ):
    return "<A HREF='%s/report_join?nid1=%s&type1=g&nid2=%s&type2=g&gop=-4&gep=-0.4&skip=0&SUBMIT=Submit+Query' target='_blank'>*</A>" % (PICASSO_URL, nid1, nid2)

def Build_join2way( nid1, nid2 ):
    return "<A HREF='%s/report_join?nid1=%s&type1=g&nid2=%s&type2=g&gop=-4&gep=-0.4&skip=1&min_score=10.0&method=0&SUBMIT=Submit+Query' target='_blank'>1->2</A>;" % (PICASSO_URL, nid1, nid2) +\
           "<A HREF='%s/report_join?nid1=%s&type1=g&nid2=%s&type2=g&gop=-4&gep=-0.4&skip=1&min_score=10.0&method=0&SUBMIT=Submit+Query' target='_blank'>2->1</A>" % (PICASSO_URL, nid2, nid1)


class RSDB2HTML_mergers_false( RSDB2HTML):

    def __init__( self ):
        RSDB2HTML.__init__( self )
        self.mTitle   = "mergers_false"
        self.mCaption = "mergers_false"

    def ParseHeader(self, line ):
        self.mHeading = ["Align"] + string.split( line[:-1], "\t")

    def ParseColumns( self, columns ):
        return [Build_join2way( columns[1], columns[10]),
                columns[0], Build_report_sequence(columns[1])] +\
                columns[2:7] +\
                [ Build_scop( columns[7] ) ] +\
                columns[8:10] +\
                [ Build_report_sequence(columns[10])] +\
                columns[11:16] +\
                [ Build_scop( columns[16] ) ] +\
                columns[17:]

class RSDB2HTML_mergers_true( RSDB2HTML):

    def __init__( self ):
        RSDB2HTML.__init__( self )
        self.mTitle   = "mergers_true"
        self.mCaption = "mergers_true"

    def ParseHeader(self, line ):
        self.mHeading = ["Align"] + string.split( line[:-1], "\t")

    def ParseColumns( self, columns ):
        return [Build_join2way( columns[1], columns[10]),
                columns[0], Build_report_sequence(columns[1])] +\
                columns[2:7] +\
                [ Build_scop( columns[7] ) ] +\
                columns[8:10] +\
                [ Build_report_sequence(columns[10])] +\
                columns[11:16] +\
                [ Build_scop( columns[16] ) ] +\
                columns[17:]
    
    
class RSDB2HTML_non_mergers_false( RSDB2HTML):

    def __init__( self ):
        RSDB2HTML.__init__( self )
        self.mTitle   = "non_mergers_false"
        self.mCaption = "non_mergers_false"

    def ParseHeader(self, line ):
        self.mHeading = ["Align"] + string.split( line[:-1], "\t")

    def ParseColumns( self, columns ):
        return [Build_join2way( columns[1], columns[10]),
                columns[0], Build_report_sequence(columns[1])] +\
                columns[2:7] +\
                [ Build_scop( columns[7] ) ] +\
                columns[8:10] +\
                [ Build_report_sequence(columns[10])] +\
                columns[11:16] +\
                [ Build_scop( columns[16] ) ] +\
                columns[17:]
    
    
    
            
        

    
