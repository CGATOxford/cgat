#!/usr/bin/env bash

# References
# http://kvz.io/blog/2013/11/21/bash-best-practices/
# http://jvns.ca/blog/2017/03/26/bash-quirks/

# exit when a command fails
set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
#set -o nounset

# trace what gets executed
#set -o xtrace

# log installation information
log() {
   echo "# install-CGAT-tools.sh log | `hostname` | `date` | $1 "
}

# message to display when the OS is not correct
sanity_check_os() {
   echo
   echo " Unsupported operating system: $OS"
   echo " Installation aborted."
   echo 
   echo " Supported operating systems are: "
   echo " Ubuntu 12.x"
   echo " CentOS 6.x"
   echo " Scientific Linux 6.x"
   echo
   exit 1;
} # sanity_check_os


# function to detect the Operating System
detect_os() {

if [[ -f /etc/os-release ]]; then

   OS=$(cat /etc/os-release | awk '/VERSION_ID/ {sub("="," "); print $2;}' | sed 's/\"//g' | awk '{sub("\\."," "); print $1;}')
   if [[ "$OS" != "12" ]] ; then

      echo       
      echo " Sorry, this version of Ubuntu has not been tested. Only Ubuntu 12.x is supported so far. "
      echo
      exit 1;

   fi

   OS="ubuntu"

elif [[ -f /etc/system-release ]]; then

   OP=$(cat /etc/system-release | awk ' {print $1;}')
   if [[ "$OP" == "Scientific" ]] ; then
      OP=$(cat /etc/system-release | awk ' {print $4;}' | awk '{sub("\\."," "); print $1;}')
      if [[ "$OP" != "6" ]] ; then
         echo
         echo " Sorry, this version of Scientific Linux has not been tested. Only 6.x versions are supported so far. "
         echo
         exit 1;
      else
         OS="sl"
      fi
   elif [[ "$OP" == "CentOS" ]] ; then
      OP=$(cat /etc/system-release | awk ' {print $3;}' | awk '{sub("\\."," "); print $1;}')
      if [[ "$OP" != "6" ]] ; then
         echo
         echo " Sorry, this version of CentOS has not been tested. Only 6.x versions are supported so far. "
         echo
         exit 1;
      else
         OS="centos"
      fi
   else
      sanity_check_os
   fi

else

   sanity_check_os

fi
} # detect_os


# install operating system dependencies
install_os_packages() {

detect_os

if [[ "$OS" == "ubuntu" ]] || [[ "$OS" == "travis" ]] ; then

   log "installing packages for Ubuntu "

   sudo apt-get --quiet install -y gcc g++ zlib1g-dev libssl-dev libssl1.0.0 libbz2-dev libfreetype6-dev libpng12-dev libblas-dev libatlas-dev liblapack-dev gfortran libpq-dev r-base-dev libreadline-dev libmysqlclient-dev libboost-dev libsqlite3-dev;

elif [[ "$OS" == "sl" ]] || [[ "$OS" == "centos" ]] ; then

   log "installing packages for Scientific Linux / CentOS "

   yum -y install gcc zlib-devel openssl-devel bzip2-devel gcc-c++ freetype-devel libpng-devel blas atlas lapack gcc-gfortran postgresql-devel R-core-devel readline-devel mysql-devel boost-devel sqlite-devel

   # additional configuration for scipy (Scientific Linux only)
   if [[ "$OS" == "sl" ]] ; then
      ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so
   fi

   # additional configuration for blas and lapack
   ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
   ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so;

else

   sanity_check_os $OS

fi # if-OS
} # install_os_packages


# detect CGAT installation
detect_cgat_installation() {

if [[ -z "$CGAT_HOME" ]] ; then

   if [[ -d "$HOME/cgat-install/conda-install" ]] ; then
      UNINSTALL_DIR="$HOME/cgat-install"
   fi

else

   if [[ -d "$CGAT_HOME/conda-install" ]] ; then
      UNINSTALL_DIR="$CGAT_HOME"
   fi

fi

} # detect_cgat_installation


# configure environment variables 
# set: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE, INSTALL_PYTHON_VERSION
get_cgat_env() {

if [[ $TRAVIS_INSTALL ]] ; then

   CGAT_HOME=$TRAVIS_BUILD_DIR
   CONDA_INSTALL_TYPE="cgat-scripts-nosetests"
   INSTALL_PYTHON_VERSION=$TRAVIS_PYTHON_VERSION

elif [[ $JENKINS_INSTALL ]] ; then

   CGAT_HOME=$WORKSPACE
   CONDA_INSTALL_TYPE="cgat-scripts-nosetests"
   INSTALL_PYTHON_VERSION=$JENKINS_PYTHON_VERSION

else

   if [[ -z $CGAT_HOME  ]] ; then
      CGAT_HOME=$HOME/cgat-install
   fi

   if [[ -z $INSTALL_PYTHON_VERSION ]] ; then
      INSTALL_PYTHON_VERSION=2
   fi

   if [[ $INSTALL_SCRIPTS ]] ; then
      CONDA_INSTALL_TYPE="cgat-scripts"
   elif [[ $INSTALL_DEVEL ]] ; then
      CONDA_INSTALL_TYPE="cgat-scripts-devel"
   elif [[ $INSTALL_TEST ]] || [[ $INSTALL_UPDATE ]] ; then
      if [[ -d $CGAT_HOME/conda-install ]] ; then
         AUX=`find $CGAT_HOME/conda-install/envs/cgat-* -maxdepth 0`
	 CONDA_INSTALL_TYPE=`basename $AUX`
      else
         echo
         echo " The location of the CGAT code was not found (function: get_cgat_env). "
	 echo " Please install it first or use --location option with full path to your installation. "
         echo
         exit 1
      fi
   else
      echo
      echo " Wrong installation type! "
      echo " Installation aborted. "
      echo
      exit 1
   fi # if install type

fi # if travis install

CONDA_INSTALL_DIR=$CGAT_HOME/conda-install
CONDA_INSTALL_ENV=$(echo $CONDA_INSTALL_TYPE | cut -c1-6)

} # get_cgat_env


# setup environment variables
setup_env_vars() {

export CFLAGS=$CFLAGS" -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib"
export CPATH=$CPATH" -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib"
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/include
export LIBRARY_PATH=$LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_ENV/lib/R/lib

} # setup_env_vars

# print related environment variables
print_env_vars() {

echo
echo " Debugging: "
echo " CFLAGS: "$CFLAGS
echo " CPATH: "$CPATH
echo " C_INCLUDE_PATH: "$C_INCLUDE_PATH
echo " CPLUS_INCLUDE_PATH: "$CPLUS_INCLUDE_PATH
echo " LIBRARY_PATH: "$LIBRARY_PATH
echo " LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
echo " CGAT_HOME: "$CGAT_HOME
echo " CONDA_INSTALL_DIR: "$CONDA_INSTALL_DIR
echo " CONDA_INSTALL_TYPE: "$CONDA_INSTALL_TYPE
echo " CONDA_INSTALL_ENV: "$CONDA_INSTALL_ENV
echo " PYTHONPATH: "$PYTHONPATH
[[ ! $INSTALL_TEST ]] && echo " INSTALL_PYTHON_VERSION: "$INSTALL_PYTHON_VERSION
[[ ! $INSTALL_TEST ]] && echo " INSTALL_BRANCH: "$INSTALL_BRANCH
echo

} # print_env_vars

# Travis installations are running out of RAM
# with large conda installations. Issue has been submitted here:
# https://github.com/conda/conda/issues/1197
# While we wait for a response, we'll try to clean up the conda
# installation folder as much as possible
conda_cleanup() {
conda clean --index-cache
conda clean --lock
conda clean --tarballs -y
conda clean --packages -y
}


# proceed with conda installation
conda_install() {

log "installing conda"

detect_cgat_installation

if [[ -n "$UNINSTALL_DIR" ]] ; then

   echo
   echo " An installation of the CGAT code was found in: $UNINSTALL_DIR"
   echo " Please use --location to install CGAT code in a different location "
   echo " or uninstall the current version before proceeding."
   echo
   echo " Installation is aborted."
   echo
   exit 1

fi

# get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env

mkdir -p $CGAT_HOME
cd $CGAT_HOME

log "downloading miniconda"
# download and install conda
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

log "installing miniconda"
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_INSTALL_DIR
export PATH="$CONDA_INSTALL_DIR/bin:$PATH"
hash -r

# install cgat environment
log "updating conda environment"
conda config --set allow_softlinks False
conda update -q conda --yes
conda info -a

log "installing conda CGAT environment"
# SLV workaround until problems are resolved for:
# * Latest R 3.3.2 available in conda: Problems with icu here https://travis-ci.org/CGATOxford/cgat/builds/240711411 and here https://github.com/ContinuumIO/anaconda-issues/issues/1403
# * Latest pysam works with the CGAT code: the gff3 support is broken, best not to use 0.11.2 for now (AH comment)
conda create -q -n $CONDA_INSTALL_ENV $CONDA_INSTALL_TYPE python=$INSTALL_PYTHON_VERSION pysam=0.11.1 r=3.3.1 gcc --override-channels --channel bioconda --channel r --channel defaults --channel conda-forge --yes

log "installing CGAT code into conda environment"
# if installation is 'devel' (outside of travis), checkout latest version from github
if [[ "$OS" != "travis" ]] ; then

   DEV_RESULT=0

   if [[ $INSTALL_DEVEL ]] ; then

      if [[ $INSTALL_ZIP ]] ; then
	 # get the latest version from Git Hub in zip format
	 cd $CGAT_HOME
         wget https://github.com/CGATOxford/cgat/archive/$INSTALL_BRANCH.zip
         unzip master.zip
	 cd cgat-master/
      else
         # get latest version from Git Hub with git clone
         git clone --branch=$INSTALL_BRANCH https://github.com/CGATOxford/cgat.git $CGAT_HOME/cgat-code
         cd $CGAT_HOME/cgat-code
      fi

      # activate cgat environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

      # SLV: workaround until bx-python is available with Python 3
      pip install bx-python

      # Set up other environment variables
      setup_env_vars

      # brute force: modify console_scripts variable/entry point for cgat command
      sed -i 's/CGATScripts/scripts/g' setup.py

      # Python preparation
      # remove install_requires (no longer required with conda package)
      sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
      sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
      python setup.py develop

      DEV_RESULT=$?

      if [[ $DEV_RESULT -ne 0 ]] ; then
         echo
         echo " There was a problem doing: 'python setup.py develop' "
         echo " Installation did not finish properly. "
         echo 
         echo " Please submit this issue via Git Hub: "
         echo " https://github.com/CGATOxford/cgat/issues "
	 echo

         print_env_vars

      fi # if-$?
   fi # if INSTALL_DEVEL

   # check whether conda create went fine
   if [[ $DEV_RESULT -ne 0 ]] ; then
      echo
      echo " There was a problem installing the code with conda. "
      echo " Installation did not finish properly. "
      echo
      echo " Please submit this issue via Git Hub: "
      echo " https://github.com/CGATOxford/cgat/issues "
      echo

      print_env_vars

   else
      clear
      echo 
      echo " The CGAT code was successfully installed!"
      echo
      echo " To activate the CGAT environment type: "
      echo " $ source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV"
      [[ $INSTALL_SCRIPTS ]] && echo " cgat --help"
      echo
      echo " To deactivate the environment, use:"
      echo " $ source deactivate"
      echo
   fi # if-$ conda create

fi # if travis install

} # conda install


# test code with conda install
conda_test() {

log "starting conda_test"

# get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env

setup_env_vars

# setup environment and run tests
if [[ $TRAVIS_INSTALL ]] || [[ $JENKINS_INSTALL ]] ; then

   # enable Conda env
   log "activating CGAT conda environment"
   source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

   # SLV: workaround until bx-python is available with Python 3
   log "pip-installing additional packages"
   pip install bx-python

   # python preparation
   log "install CGAT code into conda environment"
   cd $CGAT_HOME
   # remove install_requires (no longer required with conda package)
   sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
   sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
   python setup.py develop

   log "starting tests"
   # run nosetests
   if [[ $TEST_ALL ]] ; then
      log "test_import.py" && nosetests -v tests/test_import.py && \
      log "test_style.py" && nosetests -v tests/test_style.py && \
      echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml && \
      log "test_commandline" && nosetests -v tests/test_commandline.py && \
      log "test_scripts" && nosetests -v tests/test_scripts.py ;
   elif [[ $TEST_IMPORT ]] ; then
      nosetests -v tests/test_import.py ;
   elif [[ $TEST_STYLE ]] ; then
      nosetests -v tests/test_style.py ;
   elif [[ $TEST_CMDLINE ]] ; then
      echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml
      nosetests -v tests/test_commandline.py ;
   elif [[ $TEST_PRODUCTION_SCRIPTS  ]] ; then
      echo -e "restrict:\n    manifest:\n" > tests/_test_scripts.yaml
      nosetests -v tests/test_scripts.py ;
   else
      nosetests -v tests/test_scripts.py ;
   fi

else

   if [[ $CONDA_INSTALL_TYPE ]] ; then

      # prepare environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

      if [[ $INSTALL_ZIP ]] ; then
         cd $CGAT_HOME/cgat-master
      else
         cd $CGAT_HOME/cgat-code
      fi

      # remove install_requires (no longer required with conda package)
      sed -i'' -e '/REPO_REQUIREMENT/,/pass/d' setup.py
      sed -i'' -e '/# dependencies/,/dependency_links=dependency_links,/d' setup.py
      python setup.py develop
      OUTPUT_DIR=`pwd`

      # run tests
      /usr/bin/time -o test_import.time -v nosetests -v tests/test_import.py >& test_import.out
      if [[ $? -eq 0 ]] ; then
         echo
         echo " test_import.py passed successfully! "
         echo
      else
         echo
         echo " test_import.py failed. Please see $OUTPUT_DIR/test_import.out file for detailed output. "
         echo

         print_env_vars

      fi

      /usr/bin/time -o test_scripts.time -v nosetests -v tests/test_scripts.py >& test_scripts.out
      if [[ $? -eq 0 ]] ; then
         echo
         echo " test_scripts.py passed successfully! "
         echo
      else
         echo
         echo " test_scripts.py failed. Please see $OUTPUT_DIR/test_scripts.out file for detailed output. "
         echo

         print_env_vars

      fi
     
   else
      echo
      echo " There was an error running the tests. "
      echo " Execution aborted. "
      echo

      print_env_vars

      exit 1
   fi

fi # if-OS

} # conda_test


# update conda installation
conda_update() {

# get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env

source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV
conda update --all

if [[ ! $? -eq 0 ]] ; then

   echo
   echo " There was a problem updating the installation. "
   echo 
   echo " Please submit this issue via Git Hub: "
   echo " https://github.com/CGATOxford/cgat/issues "
   echo 

else 

   echo
   echo " All packages were succesfully updated. "
   echo 

fi 

} # conda_update


# unistall CGAT code collection
uninstall() {

detect_cgat_installation

if [[ -z "$UNINSTALL_DIR" ]] ; then

   echo
   echo " The location of the CGAT code was not found. "
   echo " Please uninstall manually."
   echo
   exit 1
    
else

   rm -rf $UNINSTALL_DIR
   if [[ $? -eq 0 ]] ; then
      echo
      echo " CGAT code successfully uninstalled."
      echo 
      exit 0
   else
      echo
      echo " There was a problem uninstalling the CGAT code."
      echo " Please uninstall manually."
      echo
      exit 1
   fi
fi

}


# function to display help message
help_message() {
echo
echo " This script uses Conda to install the CGAT Code Collection:"
echo " https://www.cgat.org/downloads/public/cgat/documentation/"
echo
echo " If you only need to use the scripts published here:"
echo "   https://www.cgat.org/downloads/public/cgat/documentation/cgat.html"
echo " type:"
echo " ./install-CGAT-tools.sh --cgat-scripts [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " The default location is: $HOME/cgat-install"
echo
echo " Otherwise, if you prefer to use the latest development version of the scripts instead, type:"
echo " ./install-CGAT-tools.sh --cgat-devel [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " Both installations create a new Conda environment ready to run the CGAT code."
echo
echo " The default Python version for CGAT is 2.7 but we are moving the code to Python 3."
echo " If you want to try our code running with Python 3, please type:"
echo " ./install-CGAT-tools.sh --cgat-devel --python 3 [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " It is also possible to install/test a specific branch of the code on github:"
echo " ./install-CGAT-tools.sh --cgat-devel --python 3 --branch <name-of-branch> [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " On the other hand, if you are looking for other alternative installation options please visit:"
echo " https://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html"
echo 
echo " To test the installation:"
echo " ./install-CGAT-tools.sh --test [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " To update the Conda packages:"
echo " ./install-CGAT-tools.sh --update [--location </full/path/to/folder/without/trailing/slash>]"
echo 
echo " To uninstall the CGAT code:"
echo " ./install-CGAT-tools.sh --uninstall [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " Please submit any issues via Git Hub:"
echo " https://github.com/CGATOxford/cgat/issues"
echo
exit 1
} # help_message

# the script starts here

if [[ $# -eq 0 ]] ; then

   help_message

fi

# these variables will store the information about input parameters
OS="default"
# travis execution
TRAVIS_INSTALL=
# jenkins testing
JENKINS_INSTALL=
# install operating system's dependencies
OS_PKGS=
# conda installation type
INSTALL_SCRIPTS=
INSTALL_DEVEL=
# test current installation
INSTALL_TEST=
# update current installation
INSTALL_UPDATE=
# uninstall CGAT code
UNINSTALL=
UNINSTALL_DIR=
# where to install CGAT code
CGAT_HOME=
# instead of cloning with git, we can download zipped CGAT code
INSTALL_ZIP=
# which python version to use
INSTALL_PYTHON_VERSION=
# which github branch to use (default: master)
INSTALL_BRANCH="master"
# variable to store input parameters
INPUT_ARGS=$(getopt -n "$0" -o htj1234567:zp:b: --long "help,
                                                        travis,
                                                        jenkins,
                                                        install-os-packages,
                                                        cgat-scripts,
                                                        cgat-devel,
                                                        test,
                                                        update,
                                                        uninstall,
                                                        location:,
                                                        zip,
                                                        python:,
                                                        branch:"  -- "$@")
eval set -- "$INPUT_ARGS"

# process all the input parameters first
while [[ "$1" != "--" ]]
do

  if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]] ; then

    help_message

  elif [[ "$1" == "--travis" ]] ; then
      
      TRAVIS_INSTALL=1
      shift ;

  elif [[ "$1" == "--jenkins" ]] ; then

      JENKINS_INSTALL=1
      shift ;

  elif [[ "$1" == "--install-os-packages" ]] ; then
      
      OS_PKGS=1
      shift ;

  elif [[ "$1" == "--cgat-scripts" ]] ; then
      
      INSTALL_SCRIPTS=1
      shift ;

  elif [[ "$1" == "--cgat-devel" ]] ; then

      INSTALL_DEVEL=1
      shift ;

  elif [[ "$1" == "--test" ]] ; then

      INSTALL_TEST=1
      shift ;

  elif [[ "$1" == "--update" ]] ; then

      INSTALL_UPDATE=1
      shift ;

  elif [[ "$1" == "--uninstall" ]] ; then

      UNINSTALL=1
      shift ; 

  elif [[ "$1" == "--location" ]] ; then

      CGAT_HOME="$2"
      shift 2 ;

  elif [[ "$1" == "--zip" ]] ; then

      INSTALL_ZIP=1
      shift ;

  elif [[ "$1" == "--python" ]] ; then

      if [[ $2 -ne 2 ]] && [[ $2 -ne 3 ]] ; then
          help_message
      else
          INSTALL_PYTHON_VERSION=$2
          shift 2 ;
      fi

  elif [[ "$1" == "--branch" ]] ; then

      INSTALL_BRANCH="$2"
      shift 2 ;

  else

    help_message

  fi # if-args
  

done # while-loop

# sanity checks
if [[ $INSTALL_SCRIPTS ]] && [[ $INSTALL_DEVEL ]] ; then
   echo
   echo " Incorrect input arguments: mixing --cgat-scripts and --cgat-devel is not permitted."
   echo " Installation aborted. Please run -h option."
   echo
   exit 1
fi

# perform actions according to the input parameters processed
if [[ $TRAVIS_INSTALL ]] ; then

  OS="travis"
  conda_install
  conda_test

elif [[ $JENKINS_INSTALL  ]] ; then

  conda_install
  conda_test

else 

  if [[ $OS_PKGS ]] ; then
     install_os_packages
  fi

  if [[ $INSTALL_SCRIPTS ]] || [[ $INSTALL_DEVEL ]] ; then
     conda_install
  fi

  if [[ $INSTALL_TEST ]] ; then
     conda_test
  fi

  if [[ $INSTALL_UPDATE ]] ; then
     conda_update
  fi

  if [[ $UNINSTALL ]] ; then
     uninstall
  fi

fi # if-variables

