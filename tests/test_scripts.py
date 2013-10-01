import subprocess
import tempfile
import os
import shutil
import re
import sys
import glob
import gzip
from nose.tools import assert_equal

######################################### 
# List of tests to perform.
######################################### 
# The fields are:
# 1. Name of the script
# 2. Filename to use as stdin
# 3. Option string
# 4. List of output files to collect
# 5. List of reference files
import yaml

def check_script( test_name, script, stdin, options, outputs, references, workingdir ):
    
    tmpdir = tempfile.mkdtemp()

    stdout = os.path.join( tmpdir, 'stdout' )
    if stdin: 
        if stdin.endswith( ".gz"):
            stdin = '< <( zcat %s/%s)' % (os.path.abspath(workingdir), stdin)
        else:
            stdin = '< %s/%s' % (os.path.abspath(workingdir), stdin)
    else: stdin = ""

    if options:
        options = re.sub( "%TMP%", tmpdir, options )
        options = re.sub( "<TMP>", tmpdir, options )
        options = re.sub( "%DIR%", os.path.abspath(workingdir), options )
        options = re.sub( "<DIR>", os.path.abspath(workingdir), options )
    else:
        options = ""

    # use /bin/bash in order to enable "<( )" syntax in shells 
    statement = ( "/bin/bash -c 'python %(script)s "
                  " %(options)s"
                  " %(stdin)s"
                  " > %(stdout)s'" ) % locals()

    retval = subprocess.call( statement, 
                              shell = True,
                              cwd = tmpdir )
    assert retval == 0

    # for version tests, do not compare output
    if test_name == "version":
        return
    
    # compare line by line, ignoring comments
    for output, reference in zip( outputs, references ):
        if output == "stdout": 
            output = stdout
        else:
            output = os.path.join( tmpdir, output )

        if not os.path.exists( output ): 
            raise OSError( "output file '%s'  does not exist" % output )
                
        reference = os.path.join( workingdir, reference)
        if not os.path.exists( reference):
            raise OSError( "reference file '%s'  does not exist (%s)" % (reference,tmpdir) )

        for a,b in zip( _read( output ), _read( reference )):
            assert_equal( a, b, "files %s and %s are not the same" % (output, reference) )

    shutil.rmtree( tmpdir )

def test_scripts():
    '''yield list of scripts to test.'''

    scriptdirs = glob.glob( "tests/*.py" )

    if os.path.exists( "tests/test_scripts.yaml" ):
        config = yaml.load( open( "tests/test_scripts.yaml" ) )
        if config != None:
            if "restrict" in config and config["restrict"]:
                values = config["restrict"]
                if "glob" in values:
                    scriptdirs = glob.glob( "tests/*.py" )
                
                if "manifest" in values:
                    # take scripts defined in the MANIFEST.in file
                    scriptdirs = [x for x in open( "MANIFEST.in" ) if x.startswith("include scripts") and x.endswith(".py\n") ]
                    scriptdirs = [re.sub( "include\s*scripts/", "tests/", x[:-1] ) for x in scriptdirs ]
                #sys.stderr.write("%s\n" % scriptdirs )

                if "regex" in values:
                    rx = re.compile( values["regex"] )
                    scriptdirs = filter( rx.search, scriptdirs )
                    

    scriptdirs.sort()        

    for scriptdir in scriptdirs:
        fn = '%s/tests.yaml' % scriptdir
        if not os.path.exists( fn ): continue

        script_tests = yaml.load( open( fn ))
        
        for test, values in script_tests.items():
            check_script.description = os.path.join( scriptdir, test)
            script_name = os.path.basename(scriptdir)

            # deal with scripts in subdirectories. These are prefixed by a "<subdir>_" 
            # for example: optic_compare_projects.py is optic/compare_procjets.py
            if "_" in script_name:
                parts = script_name.split("_")
                if os.path.exists( os.path.join( "scripts", parts[0], "_".join(parts[1:]))):
                    script_name = os.path.join(parts[0], "_".join(parts[1:]))
                    
            yield( check_script,
                   test,
                   os.path.abspath( os.path.join( "scripts", script_name )),
                   values.get('stdin', None), 
                   values['options'], 
                   values['outputs'], 
                   values['references'],
                   scriptdir )

def _read( fn ):
    if fn.endswith( ".gz" ):
        with gzip.open(fn) as inf:
            for line in inf:
                if not line.startswith("#"): yield line
    else:
        with open(fn) as inf:
            for line in inf:
                if not line.startswith("#"): yield line




