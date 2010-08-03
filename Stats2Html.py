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
Stats2Html.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
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
    
    
    
            
        

    
