"""collect scripts from a directory and create documentation
stubs.

Files starting with lower case characters are put into
the directory :file:`scripts`, while those starting with
upper case characters are put into :file:`modules`.
"""

TEMPLATE_SCRIPT='''
.. automodule:: %(prefix)s

.. program-output:: python ../scripts/%(prefix)s.py --help

'''

TEMPLATE_MODULE='''
.. automodule:: %(prefix)s
   :members:
   :show-inheritance:
'''

TEMPLATE_PIPELINE='''
.. automodule:: %(prefix)s
   :members:
   :show-inheritance:
'''

TEMPLATE_PIPELINEMODULE='''
.. automodule:: %(prefix)s
   :members:
   :show-inheritance:
'''

import glob
import os

import CGAT.Experiment as E

if __name__ == "__main__":

    E.Start()

    dirs = ( ("../scripts/*.py", TEMPLATE_SCRIPT, 'scripts'),
             ("../CGAT/*.py", TEMPLATE_MODULE, 'modules'),
             ("../CGATPipelines/pipeline*.py", TEMPLATE_PIPELINE, 'pipelines'),
             ("../CGATPipelines/[A-Z]*.py", TEMPLATE_PIPELINEMODULE, 'pipelinemodules' ) )

    ncreated, nskipped = 0, 0

    for glob_expression, template, dest in dirs:

        if not os.path.exists( dest ):
            os.mkdir( dest )

        files = glob.glob( os.path.abspath( glob_expression ) )

        for filename in files:
            dirname, name = os.path.split( filename )
            prefix = name[:-3]

            #if os.path.exists( os.path.join( dirname, "_%s.pyx" % prefix )):
            #    E.warn( "ignoring pyximport file _%s.pyx" % prefix )
            #    continue

            filename = os.path.join( os.path.abspath(dest), "%s.rst" % prefix )
            if os.path.exists( filename ):
                nskipped += 1
                continue

            E.debug( "adding %s" % filename )
            outfile = open( filename, "w" )
            outfile.write( template % locals() )
            outfile.close()

            ncreated += 1

    E.info( "ncreated=%i, nskipped=%i" % (ncreated, nskipped ) )
    
    E.Stop()
