import glob
import sys
import os
import subprocess
import re

########################################################################
#######################################################################
# Check for dependencies
#
# Is there a way to do this more elegantly?
#     1. Run "pip install numpy"
#     2. Wrap inside functions (works for numpy/pysam, but not cython)
try:
    import numpy
except ImportError:
    raise ImportError(
        "the CGAT code collection requires numpy to be installed "
        "before running setup.py (pip install numpy)")

try:
    import Cython
except ImportError:
    raise ImportError(
        "the CGAT code collection requires cython to "
        "be installed before running setup.py (pip install cython)")

try:
    import pysam
except ImportError:
    raise ImportError(
        "the CGAT code collection requires pysam to "
        "be installed before running setup.py (pip install pysam)")

########################################################################
########################################################################
# Import setuptools
# Use existing setuptools, otherwise try ez_setup.
try:
    import setuptools
except ImportError:
    # try to get via ez_setup
    # ez_setup did not work on all machines tested as
    # it uses curl with https protocol, which is not
    # enabled in ScientificLinux
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")

from Cython.Distutils import build_ext

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect CGAT version
sys.path.insert(0, "CGAT")
import version

version = version.__version__

###############################################################
###############################################################
# Check for external dependencies
#
# Not exhaustive, simply execute a representative tool from a toolkit.
external_dependencies = (
    ("wigToBigWig", "UCSC tools", 255),
    ("bedtools", "bedtools", 0),
    )

for tool, toolkit, expected in external_dependencies:
    try:
        # py3k
        from subprocess import DEVNULL
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    try:
        retcode = subprocess.call(tool, shell=True,
                                  stdout=DEVNULL, stderr=DEVNULL)
    except OSError as msg:
        print(("WARNING: depency check for %s failed: %s" % (toolkit, msg)))

    # UCSC tools return 255 when called without arguments
    if retcode != expected:
        print(("WARNING: depency check for %s(%s) failed, error %i" %
               (toolkit, tool, retcode)))

###############################################################
###############################################################
# Define dependencies
#
major, minor1, minor2, s, tmp = sys.version_info

if (major == 2 and minor1 < 7) or major < 2:
    raise SystemExit("""CGAT requires Python 2.7 or later.""")


#####################################################################
#####################################################################
# Code to install dependencies from a repository
#####################################################################
# Modified from http://stackoverflow.com/a/9125399
#####################################################################
def which(program):
    """
    Detect whether or not a program is installed.
    Thanks to http://stackoverflow.com/a/377028/70191
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

REPO_REQUIREMENT = re.compile(
    r'^-e (?P<link>(?P<vcs>git|svn|hg|bzr).+#egg=(?P<package>.+)-(?P<version>\d(?:\.\d)*))$')
HTTPS_REQUIREMENT = re.compile(
    r'^-e (?P<link>.*).+#(?P<package>.+)-(?P<version>\d(?:\.\d)*)$')
install_requires = []
dependency_links = []

for requirement in (
        l.strip() for l in open('requires.txt') if not l.startswith("#")):
    match = REPO_REQUIREMENT.match(requirement)
    if match:
        assert which(match.group('vcs')) is not None, \
            ("VCS '%(vcs)s' must be installed in order to "
             "install %(link)s" % match.groupdict())
        install_requires.append("%(package)s==%(version)s" % match.groupdict())
        dependency_links.append(match.group('link'))
        continue

    if requirement.startswith("https"):
        install_requires.append(requirement)
        continue

    match = HTTPS_REQUIREMENT.match(requirement)
    if match:
        install_requires.append("%(package)s>=%(version)s" % match.groupdict())
        dependency_links.append(match.group('link'))
        continue

    install_requires.append(requirement)

if major == 2:
    install_requires.extend(['web.py>=0.37',
                             'xlwt>=0.7.4',
                             'matplotlib-venn>=0.5'])
elif major == 3:
    pass

cgat_packages = find_packages()
cgat_package_dirs = {'CGAT': 'CGAT'}

##########################################################
##########################################################
# Classifiers
classifiers = """
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
# Cython Extensions
# Connected components cython extension
Components = Extension(
    'CGAT.Components',
    ['CGAT/Components/Components.pyx',
     'CGAT/Components/connected_components.cpp', ],
    library_dirs=[],
    libraries=[],
    language="c++",
)

# Nested containment lists
NCL = Extension(
    "CGAT.NCL.cnestedlist",
    ["CGAT/NCL/cnestedlist.pyx",
     "CGAT/NCL/intervaldb.c"],
    library_dirs=[],
    libraries=[],
    language="c",
)

# Timeseries analysis
Timeseries = Extension(
    "CGAT.Timeseries.cmetrics",
    ["CGAT/Timeseries/cmetrics.pyx"],
    include_dirs=[numpy.get_include()],
    library_dirs=[],
    libraries=[],
    language="c",
)

# Nested containment lists
GeneModelAnalysis = Extension(
    "CGAT.GeneModelAnalysis",
    ["CGAT/GeneModelAnalysis.pyx"],
    include_dirs=pysam.get_include() + [numpy.get_include()],
    library_dirs=[],
    libraries=[],
    define_macros=pysam.get_defines(),
    language="c",
)

# automatically build pyximport script extensions
pyx_files = glob.glob("CGAT/scripts/*.pyx")
script_extensions = []
pysam_dirname = os.path.dirname(pysam.__file__)
include_dirs = [numpy.get_include()] + pysam.get_include()

if IS_OSX:
    # linking against bundles does no work (and apparently is not needed)
    # within OS X
    extra_link_args = []
else:
    extra_link_args = [os.path.join(pysam_dirname, x) for x in
                       pysam.get_libraries()]

for pyx_file in pyx_files:
    script_name = os.path.basename(pyx_file)
    script_prefix = script_name[:-4]
    script_extensions.append(
        Extension("CGAT.scripts.%s" % (script_prefix),
                  sources=[pyx_file],
                  extra_link_args=extra_link_args,
                  include_dirs=include_dirs,
                  define_macros=pysam.get_defines())
    )


ext_modules = [Components, NCL, Timeseries, GeneModelAnalysis] + script_extensions

setup(
    # package information
    name='CGAT',
    version=version,
    description='CGAT : the Computational Genomics Analysis Toolkit',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='CGAT : the Computational Genomics Analysis Toolkit',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="http://www.cgat.org/cgat/Tools/",
    # package contents
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    entry_points={
        'console_scripts': ['cgat = CGAT.cgat:main']
    },
    # dependencies
    install_requires=install_requires,
    dependency_links=dependency_links,
    # extension modules
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    # other options
    zip_safe=False,
    test_suite="tests",
)
