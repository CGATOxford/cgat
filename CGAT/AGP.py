################################################################################
#   Gene prediction pipeline 
#
#   $Id: AGP.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2005 Andreas Heger
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
"""
AGP.py - working with AGP files 
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

to assemble contigs to scaffolds.

Code
----

"""

class ObjectPosition:
    def __init__(self):
        pass

    def map( self, start, end ):
        if self.mOrientation:
            return start + self.start, end + self.start
        else:
            return end + self.start, start + self.start
        
class ComponentPosition:
    def __init__(self):    
        pass


class AGP:
    
    ##------------------------------------------------------------
    def readFromFile( self, infile ):
        """read an agp file.

        Example line:

        scaffold_1      1       1199    1       W       contig_13       1       1199    +

        converts coordinates to zero-based coordinates using open/closed notation.
        In AGP nomenclature
        (http://www.ncbi.nlm.nih.gov/genome/guide/Assembly/AGP_Specification.html)
        objects (obj) like scaffolds are assembled from components (com) like contigs.

        Component types are:
        W: WGS sequence
        N: gap of specified length
        """

        self.mMapComponent2Object = {}
        self.mMapObject2Component = {}

        for line in infile:
            if line[0] == "#": continue

            data = line[:-1].split("\t")

            obj_id, obj_start, obj_end, ncoms, com_type, com_id = data[:6]

            if com_type == "N": continue                
            com_start, com_end, orientation = data[6:9]

            obj_start, obj_end = int(obj_start)-1, int(obj_end)
            com_start, com_end = int(com_start)-1, int(com_end)        

            orientation = orientation in ("+", "0", "na" )
            
            if com_start != 0:
                raise "beware, non zero com_start"

            object = ObjectPosition()
            object.mId = obj_id
            object.start = obj_start
            object.end = obj_end
            object.mOrientation = orientation
            
            self.mMapComponent2Object[com_id] = object 

    def mapLocation( self, id, start, end ):
        """map a genomic location."""

        if id not in self.mMapComponent2Object:
            raise KeyError, "id %s is not known" % (id)

        pos = self.mMapComponent2Object[id]
        return (pos.mId, ) + pos.map( start, end )

        

        
