from distribute_setup import use_setuptools
use_setuptools()

from setuptools import Extension, setup, find_packages

import glob, sys, os

major, minor1, minor2, s, tmp = sys.version_info

if major==2:
    extra_dependencies = [ 'web.py>=0.37',
                           'xlwt>=0.7.4', 
                           'matplotlib-venn>=0.5' ]
elif major==3:
    extra_dependencies = []

# Dependencies shared between python 2 and 3
shared_dependencies = [
    'sphinx>=1.0.5',
    'SphinxReport>=2.0',
    'rpy2>=2.3.4',
    'numpy>=1.7',
    'scipy>=0.11',
    'matplotlib>=1.2.1', 
    'sqlalchemy>=0.7.0', 
    'openpyxl>=1.5.7' ]

if major==2 and minor1<5 or major<2:
    raise SystemExit("""SphinxReport requires Python 2.5 or later.""")

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

setup(name='CGAT',
      version='0.1',
      description='CGAT : the CGAT code collection',
      author='Andreas Heger',
      author_email='andreas.heger@gmail.com',
      packages=find_packages(), 
      package_dir = { 'CGAT': 'CGAT',
                      'Pipelines' : 'CGATPipelines' },
      url="http://code.google.com/p/sphinx-report/",
      package_data={'glob.glob': ['./templates/*', './images/*']},
      scripts=glob.glob( 'scripts/*.py' ) +\
          glob.glob( 'CGATPipelines/pipeline*.py'),
      license="BSD",
      platforms=["any",],
      keywords="report generator sphinx matplotlib sql",
      long_description='CGAT : the CGAT code collection',
      classifiers = filter(None, classifiers.split("\n")),
      install_requires = shared_dependencies + extra_dependencies,
      zip_safe = False,
      include_package_data = True,
      )

