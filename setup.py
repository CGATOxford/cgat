import glob
import sys
import os

#######################################################################
## Check for dependencies
##
## Is there a way to do this more elegantly?
##     1. Run "pip install numpy"
##     2. Wrap inside functions (works for numpy/pysam, but not cython)
#######################################################################
#######################################################################
#######################################################################
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
## Use existing setuptools
try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    ## try to get via ez_setup
    ## ez_setup did not work on all machines tested as
    ## it uses curl with https protocol, which is not
    ## enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import setup, find_packages, Extension

from Cython.Distutils import build_ext

###############################################################
###############################################################
# Perform a CGAT Code Collection Installation
INSTALL_CGAT = True

major, minor1, minor2, s, tmp = sys.version_info

if major==2:
    extra_dependencies = [ 'web.py>=0.37',
                           'xlwt>=0.7.4', 
                           'matplotlib-venn>=0.5' ]
elif major==3:
    extra_dependencies = []

if major==2 and minor1<7 or major<2:
    raise SystemExit("""CGAT requires Python 2.6 or later.""")

# Dependencies shared between python 2 and 3
shared_dependencies = [
    'sphinx>=1.0.5',
    'SphinxReport>=2.0',
    'rpy2>=2.3.4',
    'numpy>=1.7',
    'scipy>=0.11',
    'matplotlib>=1.2.1', 
    'sqlalchemy>=0.7.0', 
    'pysam>=0.7',
    'openpyxl>=1.5.7',
    'MySQL-python>1.2.3',
    'biopython>=1.61',
    'scipy>=0.7.0',
    # Out of date, install manually latest version
    # 'bx-python>=0.7.1',
    'networkx>=1.8.1',
    'PyGreSQL>=4.1.1',
    'drmaa>=0.5',
    'ruffus>=2.2',
    'pybedtools>=0.6.2',
    'rdflib>=0.4.1',
    'hgapi>=1.3.0',
    'threadpool>=1.2.7',
    'PyYAML>=3.1.0',
    'pandas>=0.12.0',
    'weblogo>=3.0',
    'sphinxcontrib-programoutput>=0.8',
    'alignlib>=0.1']

cgat_packages= find_packages( exclude=["CGATPipelines*", "scripts*"])
# rename scripts to CGATScripts
cgat_packages.append( "CGATScripts" )

cgat_package_dirs = { 'CGAT': 'CGAT',
                      'CGATScripts' : 'scripts' }

classifiers="""
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

# external dependencies
# R
# sqlite
# R - ggplot2
# R - RSqlite
# R - gplots (for r-heatmap)
# graphvis - for dependency graphs in documentation

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
for pyx_file in pyx_files:
    script_name = os.path.basename( pyx_file )
    script_prefix = script_name[:-4]
    script_extensions.append( 
        Extension( "CGAT.%s" % (script_prefix),
                   sources = [pyx_file],
                   extra_link_args=[ os.path.join( pysam_dirname, "csamtools.so")],
                   include_dirs = include_dirs,
                   define_macros = pysam.get_defines() )
        )

setup(## package information
    name='CGAT',
    version='0.1.4',
    description='CGAT : the Computational Genomics Analysis Toolkit',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="BSD",
    platforms=["any",],
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
    install_requires=shared_dependencies + extra_dependencies, 
    ## extension modules
    ext_modules=[Components, NCL, Nubiscan] + script_extensions,
    cmdclass = {'build_ext': build_ext},
    ## other options
    zip_safe = False,
    test_suite = "tests",
    )

