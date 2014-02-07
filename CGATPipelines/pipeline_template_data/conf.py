##############################################
##
## conf.py - configuration options for sphinx
##           
## This file execfile's a generic conf.py
## in the CGAT code collection.
##
#############################################
import os
import CGATPipelines

conf_file = os.path.join( os.path.dirname(CGATPipelines.__file__),
                          'configuration',
                          'conf.py' )
execfile( conf_file )
