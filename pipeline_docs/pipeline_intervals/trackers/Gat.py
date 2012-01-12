import os, sys, re, types, itertools, math

from IntervalReport import *

class GatTracker( IntervalTracker ):
    pass

class GatGenomicContextTable( GatTracker ):
    pattern = "gat_context_(.*)$"

    def __call__(self, track):
        return self.getAll( "SELECT * FROM gat_context_%(track)s" )

class GatGenomicAnnotationTable( GatTracker ):
    pattern = "gat_annotation_(.*)$"

    def __call__(self, track):
        return self.getAll( "SELECT * FROM gat_annotation_%(track)s" )
