import os
import sys
import re
import types
from VariantsReport import *


class PantherTracker(VariantsTracker):

    mPattern = "^panther$"
    mAsTables = True


class PantherCounts(StrainTracker):

    '''returns overview of panther run.
    '''

    tablename = "panther"
    tablename_map = "polyphen_map"

    def getSlices(self, subset=None):
        return ["snps", "proteins", "loci"]

    def __call__(self, track, slice=None):

        data = odict()
        if slice == "snps":
            s = "DISTINCT map.snp_id"
        elif slice == "proteins":
            s = "DISTINCT map.protein_id"
        elif slice == "loci":
            s = "DISTINCT map.locus_id"
        else:
            s = "*"

        if track == "all":
            where = "1"
        else:
            where = "map.track = '%(track)s'" % locals()

        data["total"] = self.getValue( '''SELECT COUNT(%(s)s) FROM %(tablename)s AS t, %(tablename_map)s as map
                                            WHERE map.snp_id = t.snp_id 
                                                  AND %(where)s
                                           ''' % self.members(locals() ))
        data.update( odict( self.get( '''SELECT message, COUNT(%(s)s) 
                                            FROM %(tablename)s AS t, %(tablename_map)s as map
                                            WHERE map.snp_id = t.snp_id 
                                                  AND %(where)s
                                            GROUP BY message''' % self.members(locals() ) ) ))

        data["predicted"] = data[None]
        del data[None]

        return data


class PantherResults(StrainTracker):

    '''returns overview of panther results.'''

    tablename = "panther"
    tablename_map = "polyphen_map"
    pdeleterious = 0.5

    def __call__(self, track, slice=None):

        data = odict()

        s = "DISTINCT locus_id"

        if track == "all":
            where = "1"
        else:
            where = "map.track = '%(track)s'" % locals()

        data["deleterious"] = self.getValue( '''
        SELECT COUNT(%(s)s) 
        FROM %(tablename)s AS t, %(tablename_map)s as map
        WHERE map.snp_id = t.snp_id AND %(where)s AND Pdeleterious > %(pdeleterious)f ''' % self.members(locals() ) )

        data["neutral"] = self.getValue( '''
        SELECT COUNT(%(s)s) 
        FROM %(tablename)s AS t, %(tablename_map)s as map
        WHERE map.snp_id = t.snp_id AND %(where)s AND Pdeleterious <= %(pdeleterious)f ''' % self.members(locals() ) )

        return data


class PantherDistribution(PantherTracker):

    '''return distributions for various columns.'''

    def getSlices(self, subset=None):
        return ("Pdeleterious", "Nic", "subPSEC")

    def __call__(self, track, slice=None):

        return odict((
            (slice, self.getValues("SELECT %(slice)s FROM %(track)s WHERE message IS NULL" % locals())), ))
