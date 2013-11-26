import glob
import sys
import os
import subprocess

########################################################################
#######################################################################
## Check for dependencies
##
## Is there a way to do this more elegantly?
##     1. Run "pip install numpy"
##     2. Wrap inside functions (works for numpy/pysam, but not cython)
try:
    import numpy
except ImportError:
    raise ImportError("the CGAT code collection requires numpy to be installed before running setup.py (pip install numpy)" )

try:
    import Cython
except ImportError:
    raise ImportError("the CGAT code collection requires cython to be installed before running setup.py (pip install cython)" )

try: 
    import pysam
except ImportError:
    raise ImportError("the CGAT code collection requires pysam to be installed before running setup.py (pip install pysam)" )

########################################################################
########################################################################
## Import setuptools
## Use existing setuptools, otherwise try ez_setup.
try:
    import setuptools
except ImportError:
    ## try to get via ez_setup
    ## ez_setup did not work on all machines tested as
    ## it uses curl with https protocol, which is not
    ## enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion, StrictVersion
if LooseVersion( setuptools.__version__ ) < LooseVersion( '1.1' ):
    raise ImportError("the CGAT code collection requires setuptools 1.1 higher")

from Cython.Distutils import build_ext

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect CGAT version
sys.path.insert( 0, "scripts")
import version

version = version.__version__

###############################################################
###############################################################
# Check for external dependencies
#
# Not exhaustive, simply execute a representative tool from a toolkit.
external_dependencies = ( 
    ("wigToBigWig", "UCSC tools", 255),
    ("bedtools", "bedtools", 0 ),
    )

for tool, toolkit, expected in external_dependencies:
    try:
        from subprocess import DEVNULL # py3k
    except ImportError:
        import os
        DEVNULL = open(os.devnull, 'wb')

    try:
        retcode = subprocess.call( tool, shell = True, stdout=DEVNULL, stderr=DEVNULL )
    except OSError, msg:
        print( "WARNING: depency check for %s failed: %s" % (toolkit, msg ) )

    # UCSC tools return 255 when called without arguments
    if retcode != expected:
        print ( "WARNING: depency check for %s(%s) failed, error %i" % \
                    (toolkit, tool, retcode ) )

###############################################################
###############################################################
# Define dependencies 
#
# Perform a CGAT Code Collection Installation
INSTALL_CGAT_CODE_COLLECTION = True

major, minor1, minor2, s, tmp = sys.version_info

if major==3:
    raise SystemExit("""CGAT is not fully python3 compatible""")

if (major==2 and minor1<7) or major<2:
    raise SystemExit("""CGAT requires Python 2.7 or later.""")


dependencies = open( "requires.txt" ).readlines() 
dependencies = [x[:-1] for x in dependencies if x.startswith("#")]

if major==2:
    dependencies.extend( [ 'web.py>=0.37',
                           'xlwt>=0.7.4', 
                           'matplotlib-venn>=0.5' ] )
elif major==3:
    pass

if INSTALL_CGAT_CODE_COLLECTION:
    cgat_packages= find_packages( exclude=["CGATPipelines*", "scripts*"])
else:
    cgat_packages= find_packages( exclude=["scripts*"])

# rename scripts to CGATScripts
cgat_packages.append( "CGATScripts" )

cgat_package_dirs = { 'CGAT': 'CGAT',
                      'CGATScripts' : 'scripts',
                      'CGATPipelines': 'CGATPipelines' }

##########################################################
##########################################################
# Classifiers
classifiers="""
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

##########################################################
##########################################################
## Extensions
# Connected components cython extension
Components = Extension(
     'CGAT.Components',                   
     [ 'CGAT/Components/Components.pyx',        
       'CGAT/Components/connected_components.cpp',],
     library_dirs=[],
     libraries=[],              
     language="c++",               
     )

# Nested containment lists
NCL = Extension(
    "CGAT.NCL.cnestedlist",                   
    [ "CGAT/NCL/cnestedlist.pyx",    
      "CGAT/NCL/intervaldb.c" ],
      library_dirs=[],
      libraries=[],
      language="c",
    )

# Nubiscan motif mapping
Nubiscan = Extension(
    "CGAT.Nubiscan.cnubiscan",                   
    [ 'CGAT/Nubiscan/cnubiscan.pyx'],
    library_dirs=[],
    libraries=[],            
    include_dirs = [numpy.get_include()], 
    language="c",               
 )

# Automatically build script extensions
pyx_files = glob.glob( "scripts/*.pyx" )
script_extensions = []
pysam_dirname = os.path.dirname( pysam.__file__ )
include_dirs = [ numpy.get_include() ] + pysam.get_include()

if IS_OSX:
    # linking against bundles does no work (and apparently is not needed)
    # within OS X
    extra_link_args = []
else:
    extra_link_args = [ os.path.join( pysam_dirname, "csamtools.so")]

for pyx_file in pyx_files:
    script_name = os.path.basename( pyx_file )
    script_prefix = script_name[:-4]
    script_extensions.append( 
        Extension( "CGAT.%s" % (script_prefix),
                   sources = [pyx_file],
                   extra_link_args = extra_link_args,
                   include_dirs = include_dirs,
                   define_macros = pysam.get_defines() )
        )

setup(## package information
    name = 'CGAT',
    version= version,
    description = 'CGAT : the Computational Genomics Analysis Toolkit',
    author = 'Andreas Heger',
    author_email = 'andreas.heger@gmail.com',
    license = "BSD",
    platforms = ["any",],
    keywords="computational genomics",
    long_description='CGAT : the Computational Genomics Analysis Toolkit',
    classifiers = filter(None, classifiers.split("\n")),
    url="http://www.cgat.org/cgat/Tools/",
    ## package contents
    packages = cgat_packages, 

    package_dir= cgat_package_dirs,
    # package_data = { 'CGATScripts' : ['./scripts/*.py', './scripts/*.pyx', 
    #                                   './scripts/*.pyxbld', './scripts/*.pl' ],
    #                },

    include_package_data = True,

    entry_points = {
        'console_scripts': ['cgat = CGATScripts.cgat:main' ]
        },

    ## dependencies
    install_requires=dependencies,
    ## extension modules
    ext_modules=[Components, NCL, Nubiscan] + script_extensions,
    cmdclass = {'build_ext': build_ext},
    ## other options
    zip_safe = False,
    test_suite = "tests",
    )

