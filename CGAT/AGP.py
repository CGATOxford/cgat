"""AGP.py - working with AGP files 
=====================================================

This module contains a parser for reading from :term:`agp` formatted
files.

Code
----

"""


class ObjectPosition(object):

    def map(self, start, end):
        if self.mOrientation:
            return start + self.start, end + self.start
        else:
            return end + self.start, start + self.start


class AGP(object):
    """Parser for AGP formatted files."""

    def readFromFile(self, infile):
        """read an agp file.

        Example line::

           scaffold_1      1       1199    1       W       contig_13       1       1199    +

        This method converts coordinates to zero-based coordinates
        using open/closed notation.

        In AGP nomenclature
        (http://www.ncbi.nlm.nih.gov/genome/guide/Assembly/AGP_Specification.html)
        objects (obj) like scaffolds are assembled from components
        (com) like contigs.

        Component types are:

        W
           WGS sequence
        N
           gap of specified length.

        """

        self.mMapComponent2Object = {}
        self.mMapObject2Component = {}

        for line in infile:
            if line[0] == "#":
                continue

            data = line[:-1].split("\t")

            obj_id, obj_start, obj_end, ncoms, com_type, com_id = data[:6]

            if com_type == "N":
                continue
            com_start, com_end, orientation = data[6:9]

            obj_start, obj_end = int(obj_start) - 1, int(obj_end)
            com_start, com_end = int(com_start) - 1, int(com_end)

            orientation = orientation in ("+", "0", "na")

            if com_start != 0:
                raise ValueError("non zero com_start")

            object = ObjectPosition()
            object.mId = obj_id
            object.start = obj_start
            object.end = obj_end
            object.mOrientation = orientation

            self.mMapComponent2Object[com_id] = object

    def mapLocation(self, id, start, end):
        """map a genomic location.

        Raises
        ------

        KeyError
           If `id` is not present.

        """

        if id not in self.mMapComponent2Object:
            raise KeyError("id %s is not known" % (id))

        pos = self.mMapComponent2Object[id]
        return (pos.mId, ) + pos.map(start, end)
