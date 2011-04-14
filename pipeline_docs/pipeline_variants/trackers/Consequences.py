import os, sys, re, types
from VariantsReport import *

class ConsequencesTracker( VariantsTracker ):

    mPattern = "_annotation$"
    tablename_map = "polyphen_map"

class Substitutions( VariantsTracker ):
    pass
