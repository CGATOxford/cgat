#!/usr/bin/env bash

# References
# http://kvz.io/blog/2013/11/21/bash-best-practices/
# http://jvns.ca/blog/2017/03/26/bash-quirks/

# exit when a command fails
#set -o errexit

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
   CONDA_INSTALL_TYPE="scripts-nosetests.yml"
   INSTALL_PYTHON_VERSION=$TRAVIS_PYTHON_VERSION

elif [[ $JENKINS_INSTALL ]] ; then

   CGAT_HOME=$WORKSPACE
   CONDA_INSTALL_TYPE="scripts-nosetests.yml"
   INSTALL_PYTHON_VERSION=$JENKINS_PYTHON_VERSION

else

   if [[ -z $CGAT_HOME  ]] ; then
      CGAT_HOME=$HOME/cgat-install
   fi

   if [[ -z $INSTALL_PYTHON_VERSION ]] ; then
      INSTALL_PYTHON_VERSION=3
   fi

   if [[ $INSTALL_SCRIPTS ]] ; then
      CONDA_INSTALL_TYPE="cgat-scripts"
   elif [[ $INSTALL_DEVEL ]] ; then
      CONDA_INSTALL_TYPE="scripts-devel.yml"
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
CONDA_INSTALL_ENV="cgat-s"

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
source ${CONDA_INSTALL_DIR}/bin/activate
hash -r

# install cgat environment
log "updating conda environment"
conda update -q conda --yes
conda info -a

log "installing CGAT environment"
# Now using conda environment files:
# https://conda.io/docs/using/envs.html#use-environment-from-file

if [[ $INSTALL_SCRIPTS ]] ; then
   conda create -q -n ${CONDA_INSTALL_ENV} cgat-scripts --override-channels -c bioconda -c conda-forge -c defaults -y
else
   # to-update once we merge to master:
   # wget -O env.yml https://raw.githubusercontent.com/CGATOxford/cgat/conda/environments/${CONDA_INSTALL_TYPE}
   wget -O env.yml https://raw.githubusercontent.com/CGATOxford/cgat/SLV-update-install/conda/environments/${CONDA_INSTALL_TYPE}
   conda env create -f env.yml
fi

log "installing CGAT code into conda environment"
# if installation is 'devel' (outside of travis), checkout latest version from github
if [[ "$OS" != "travis" ]] ; then

   DEV_RESULT=0

   if [[ $INSTALL_DEVEL ]] ; then

      # make sure you are in the CGAT_HOME folder
      cd $CGAT_HOME

      if [[ $INSTALL_ZIP ]] ; then
	 # get the latest version from Git Hub in zip format
         wget https://github.com/CGATOxford/cgat/archive/$INSTALL_BRANCH.zip
         unzip master.zip
         rm master.zip
         mv cgat-master/ cgat-scripts/
      else
         # get latest version from Git Hub with git clone
         git clone --branch=$INSTALL_BRANCH https://github.com/CGATOxford/cgat.git $CGAT_HOME/cgat-scripts
      fi

      # make sure you are in the CGAT_HOME/cgat-scripts folder
      cd $CGAT_HOME/cgat-scripts

      # activate cgat environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV

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
      echo " The code was successfully installed!"
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

   # show conda environment used for testing
   conda env export

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

   source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_ENV
   conda list | grep cgat-scripts > /dev/null

   if [[ $? -eq 1 ]] ; then
      # this is "cgat-devel" so tests can be run

      # make sure you are in the CGAT_HOME/cgat-scripts folder
      cd $CGAT_HOME/cgat-scripts

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
      # in this case, the installation found was "cgat-scripts" so no need to run tests
      echo
      echo " You installed the cgat-scripts, which has been properly tested before. "
      echo " No need to test. Exiting now... "
      echo

      exit 0
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
echo " It is also possible to install/test a specific branch of the code on github:"
echo " ./install-CGAT-tools.sh --cgat-devel --branch <name-of-branch> [--location </full/path/to/folder/without/trailing/slash>]"
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

