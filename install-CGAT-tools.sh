#!/usr/bin/env bash

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

if [ -f /etc/os-release ]; then

   OS=$(cat /etc/os-release | awk '/VERSION_ID/ {sub("="," "); print $2;}' | sed 's/\"//g' | awk '{sub("\\."," "); print $1;}')
   if [ "$OS" != "12" ] ; then

      echo       
      echo " Sorry, this version of Ubuntu has not been tested. Only Ubuntu 12.x is supported so far. "
      echo
      exit 1;

   fi

   OS="ubuntu"

elif [ -f /etc/system-release ]; then

   OP=$(cat /etc/system-release | awk ' {print $1;}')
   if [ "$OP" == "Scientific" ] ; then
      OP=$(cat /etc/system-release | awk ' {print $4;}' | awk '{sub("\\."," "); print $1;}')
      if [ "$OP" != "6" ] ; then
         echo
         echo " Sorry, this version of Scientific Linux has not been tested. Only 6.x versions are supported so far. "
         echo
         exit 1;
      else
         OS="sl"
      fi
   elif [ "$OP" == "CentOS" ] ; then
      OP=$(cat /etc/system-release | awk ' {print $3;}' | awk '{sub("\\."," "); print $1;}')
      if [ "$OP" != "6" ] ; then
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

if [ "$OS" == "ubuntu" ] || [ "$OS" == "travis" ] ; then

   echo
   echo " Installing packages for Ubuntu "
   echo

   sudo apt-get --quiet install -y gcc g++ zlib1g-dev libssl-dev libssl1.0.0 libbz2-dev libfreetype6-dev libpng12-dev libblas-dev libatlas-dev liblapack-dev gfortran libpq-dev r-base-dev libreadline-dev libmysqlclient-dev libboost-dev libsqlite3-dev;

elif [ "$OS" == "sl" ] || [ "$OS" == "centos" ] ; then

   echo 
   echo " Installing packages for Scientific Linux / CentOS "
   echo

   yum -y install gcc zlib-devel openssl-devel bzip2-devel gcc-c++ freetype-devel libpng-devel blas atlas lapack gcc-gfortran postgresql-devel R-core-devel readline-devel mysql-devel boost-devel sqlite-devel

   # additional configuration for scipy (Scientific Linux only)
   if [ "$OS" == "sl" ] ; then
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

if [ -z "$CGAT_HOME" ] ; then

   if [ -d "$HOME/cgat-install/conda-install" ] ; then
      UNINSTALL_DIR="$HOME/cgat-install"
   fi

else

   if [ -d "$CGAT_HOME/conda-install" ] ; then
      UNINSTALL_DIR="$CGAT_HOME"
   fi

fi

} # detect_cgat_installation


# configure environment variables 
# set: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env() {

if [ $TRAVIS_INSTALL ] ; then

   CGAT_HOME=$TRAVIS_BUILD_DIR
   CONDA_INSTALL_TYPE="cgat-devel"

else

   if [ -z $CGAT_HOME  ] ; then
      CGAT_HOME=$HOME/cgat-install
   fi

   if [ "$INSTALL_SCRIPTS" == "1" ] ; then
      if [ "$INSTALL_LITE" == "1" ] ; then
         CONDA_INSTALL_TYPE="cgat-scripts-lite"
      elif [ "$INSTALL_DEVEL" == "1" ] ; then
         CONDA_INSTALL_TYPE="cgat-scripts"
      else
         CONDA_INSTALL_TYPE="cgat-scripts"
      fi
   elif [ "$INSTALL_DEVEL" == "1" ] ; then
      if [ "$INSTALL_LITE" == "1" ] ; then
         CONDA_INSTALL_TYPE="cgat-devel-lite"
      elif [ "$INSTALL_DEVEL" == "1" ] ; then
         CONDA_INSTALL_TYPE="cgat-devel"
      else
         CONDA_INSTALL_TYPE="cgat-devel"
      fi
   elif [ $INSTALL_TEST ] || [ $INSTALL_UPDATE ] ; then
      if [ -d $CGAT_HOME/conda-install ] ; then
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

} # get_cgat_env


# setup environment variables
setup_env_vars() {

export CFLAGS=$CFLAGS" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
export CPATH=$CPATH" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib/R/lib

} # setup_env_vars


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

detect_cgat_installation

if [ -n "$UNINSTALL_DIR" ] ; then

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

cd $CGAT_HOME

# download and install conda
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p $CONDA_INSTALL_DIR
export PATH="$CONDA_INSTALL_DIR/bin:$PATH"
hash -r

# install cgat environment
conda update -q conda --yes
conda info -a

# keep rpy2-2.4 for production scripts
if [ "$CONDA_INSTALL_TYPE" == "cgat-scripts" ] ; then

   conda create -q -n $CONDA_INSTALL_TYPE $CONDA_INSTALL_TYPE gcc=4.8.3 rpy2=2.4 --override-channels --channel https://conda.binstar.org/cgat --channel defaults --channel https://conda.binstar.org/r --yes

else

   conda create -q -n $CONDA_INSTALL_TYPE $CONDA_INSTALL_TYPE gcc=4.8.3 --override-channels --channel https://conda.binstar.org/cgat --channel defaults --channel https://conda.binstar.org/r --yes

fi

# if installation is 'devel' (outside of travis), checkout latest version from github
if [ "$OS" != "travis" ] ; then

   DEV_RESULT=0

   if [ $INSTALL_DEVEL ] ; then

      if [ $INSTALL_ZIP ] ; then
	 # get the latest version from Git Hub in zip format
	 cd $CGAT_HOME
         wget --no-check-certificate https://github.com/CGATOxford/cgat/archive/master.zip
         unzip master.zip
	 cd cgat-master/
      else
         # get latest version from Git Hub with git clone
         git clone https://github.com/CGATOxford/cgat.git $CGAT_HOME/cgat-code
         cd $CGAT_HOME/cgat-code
      fi

      # activate cgat environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE

      # Set up other environment variables
      setup_env_vars

      # brute force: modify console_scripts variable/entry point for cgat command
      sed -i 's/CGATScripts/scripts/g' setup.py

      # Python preparation
      python setup.py develop

      DEV_RESULT=$?

      if [ $DEV_RESULT -ne 0 ] ; then
         echo
         echo " There was a problem doing: 'python setup.py develop' "
         echo " Installation did not finish properly. "
         echo 
         echo " Please submit this issue via Git Hub: "
         echo " https://github.com/CGATOxford/cgat/issues "
	 echo
         echo " Debugging: "
         echo " CFLAGS: "$CFLAGS
         echo " CPATH: "$CPATH
         echo " C_INCLUDE_PATH: "$C_INCLUDE_PATH
         echo " LIBRARY_PATH: "$LIBRARY_PATH
         echo " LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
         echo " PYTHONPATH: "$PYTHONPATH
         echo 
      fi # if-$?
   fi # if INSTALL_DEVEL

   # check whether conda create went fine
   if [ $DEV_RESULT -ne 0 ] ; then
      echo
      echo " There was a problem installing the code with conda. "
      echo " Installation did not finish properly. "
      echo
      echo " Please submit this issue via Git Hub: "
      echo " https://github.com/CGATOxford/cgat/issues "
      echo
      echo " Debugging: "
      echo " CGAT_HOME: "$CGAT_HOME
      echo " CONDA_INSTALL_DIR: "$CONDA_INSTALL_DIR
      echo " CONDA_INSTALL_TYPE: "$CONDA_INSTALL_TYPE
      echo
   else
      clear
      echo 
      echo " The CGAT code was successfully installed!"
      echo
      echo " To activate the CGAT environment type: "
      echo " $ source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE"
      [ $INSTALL_SCRIPTS ] && echo " cgat --help"
      echo
      echo " To deactivate the environment, use:"
      echo " $ source deactivate"
      echo
   fi # if-$ conda create

fi # if travis install

} # conda install


# test code with conda install
conda_test() {

# get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env

setup_env_vars

# setup environment and run tests
if [ $TRAVIS_INSTALL ] ; then

   # enable Conda env
   source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE

   # python preparation
   cd $CGAT_HOME
   python setup.py develop

   # run nosetests
   if [ $TEST_IMPORT ] ; then
      nosetests -v tests/test_import.py ;
   elif [ $TEST_STYLE ] ; then
      nosetests -v tests/test_style.py ;
   elif [ $TEST_CMDLINE ] ; then
      echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml
      nosetests -v tests/test_commandline.py ;
   elif [ $TEST_PRODUCTION_SCRIPTS  ] ; then
      echo -e "restrict:\n    manifest:\n" > tests/_test_scripts.yaml
      nosetests -v tests/test_scripts.py ;
   else
      nosetests -v tests/test_scripts.py ;
   fi

else

   if [ "$CONDA_INSTALL_TYPE" == "cgat-scripts-lite" ] || [ "$CONDA_INSTALL_TYPE" == "cgat-scripts" ] ; then
      echo
      echo " You are using the CGAT Code Collection uploaded to pip. "
      echo " This version of the code has been well tested before release. "
      echo " Nothing to test. "
      echo
   elif [ "$CONDA_INSTALL_TYPE" == "cgat-devel-lite" ] || [ "$CONDA_INSTALL_TYPE" == "cgat-devel" ] ; then
      # prepare environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE

      if [ $INSTALL_ZIP ] ; then
         cd $CGAT_HOME/cgat-master
      else
         cd $CGAT_HOME/cgat-code
      fi

      python setup.py develop
      OUTPUT_DIR=`pwd`

      # run tests
      /usr/bin/time -o test_import.time -v nosetests -v tests/test_import.py >& test_import.out
      if [ $? -eq 0 ] ; then
         echo
         echo " test_import.py passed successfully! "
         echo
      else
         echo
         echo " test_import.py failed. Please see $OUTPUT_DIR/test_import.out file for detailed output. "
         echo
         echo " Debugging: "
         echo " CFLAGS: "$CFLAGS
         echo " CPATH: "$CPATH
         echo " C_INCLUDE_PATH: "$C_INCLUDE_PATH
         echo " LIBRARY_PATH: "$LIBRARY_PATH
         echo " LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
         echo " PYTHONPATH: "$PYTHONPATH
         echo
      fi

      /usr/bin/time -o test_scripts.time -v nosetests -v tests/test_scripts.py >& test_scripts.out
      if [ $? -eq 0 ] ; then
         echo
         echo " test_scripts.py passed successfully! "
         echo
      else
         echo
         echo " test_scripts.py failed. Please see $OUTPUT_DIR/test_scripts.out file for detailed output. "
         echo
         echo " Debugging: "
         echo " CFLAGS: "$CLAFGS
         echo " CPATH: "$CPATH
         echo " C_INCLUDE_PATH: "$C_INCLUDE_PATH
         echo " LIBRARY_PATH: "$LIBRARY_PATH
         echo " LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
         echo " PYTHONPATH: "$PYTHONPATH
         echo
      fi
     
   else
      echo
      echo " There was an error running the tests. "
      echo " Execution aborted. "
      echo
      echo " Debugging: "
      echo " CONDA_INSTALL_DIR: "$CONDA_INSTALL_DIR
      echo " CONDA_INSTALL_TYPE: "$CONDA_INSTALL_TYPE
      echo " CGAT_HOME: "$CGAT_HOME
      echo
      exit 1
   fi

fi # if-OS

} # conda_test


# update conda installation
conda_update() {

# get environment variables: CGAT_HOME, CONDA_INSTALL_DIR, CONDA_INSTALL_TYPE
get_cgat_env

source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE
conda update --all

if [ ! $? -eq 0 ] ; then

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

if [ -z "$UNINSTALL_DIR" ] ; then

   echo
   echo " The location of the CGAT code was not found. "
   echo " Please uninstall manually."
   echo
   exit 1
    
else

   rm -rf $UNINSTALL_DIR
   if [ $? -eq 0 ] ; then
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
echo " Otherwise, if you prefer to use the scripts and the pipelines altogether instead, type:"
echo " ./install-CGAT-tools.sh --cgat-devel [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " Both installations create a new Conda environment ready to run the CGAT code."
echo " On the other hand, if you are looking for other advanced installation options please visit:"
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

if [ $# -eq 0 ] ; then

   help_message

fi

# these variables will store the information about input parameters
# travis execution
TRAVIS_INSTALL=
# install operating system's dependencies
OS_PKGS=
# conda installation type
INSTALL_SCRIPTS=
INSTALL_DEVEL=
# is installation lite or full?
INSTALL_LITE=
INSTALL_FULL=
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
# variable to store input parameters
INPUT_ARGS=$(getopt -n "$0" -o h0123456789:z --long "help,
                                                  travis,
                                                  install-os-packages,
                                                  cgat-scripts,
                                                  cgat-devel,
                                                  lite,
                                                  full,
                                                  test,
                                                  update,
						  uninstall,
                                                  location:,
						  zip"  -- "$@")
eval set -- "$INPUT_ARGS"

# process all the input parameters first
while [ "$1" != "--" ]
do

  if [ "$1" == "-h" ] || [ "$1" == "--help" ] ; then

    help_message

  elif [ "$1" == "--travis" ] ; then
      
      TRAVIS_INSTALL=1
      shift ;

  elif [ "$1" == "--install-os-packages" ] ; then
      
      OS_PKGS=1
      shift ;

  elif [ "$1" == "--cgat-scripts" ] ; then
      
      INSTALL_SCRIPTS=1
      shift ;

  elif [ "$1" == "--cgat-devel" ] ; then

      INSTALL_DEVEL=1
      shift ;

  elif [ "$1" == "--lite" ] ; then

      INSTALL_LITE=1
      shift ;

  elif [ "$1" == "--full" ] ; then

      INSTALL_FULL=1
      shift ;

  elif [ "$1" == "--test" ] ; then

      INSTALL_TEST=1
      shift ;

  elif [ "$1" == "--update" ] ; then

      INSTALL_UPDATE=1
      shift ;

  elif [ "$1" == "--uninstall" ] ; then

      UNINSTALL=1
      shift ; 

  elif [ "$1" == "--location" ] ; then

      CGAT_HOME="$2"
      shift 2 ;

  elif [ "$1" == "--zip" ] ; then

      INSTALL_ZIP=1
      shift ;

  else

    help_message

  fi # if-args
  

done # while-loop

# sanity checks
if [ $INSTALL_LITE ] && [ $INSTALL_FULL ] ; then

   echo 
   echo " Incorrect input arguments: mixing --full and --lite options is not permitted."
   echo " Installation aborted. Please run -h option."
   echo
   exit 1

elif [ $INSTALL_SCRIPTS ] && [ $INSTALL_DEVEL ] ; then

   echo
   echo " Incorrect input arguments: mixing --cgat-scripts and --cgat-devel is not permitted."
   echo " Installation aborted. Please run -h option."
   echo
   exit 1

fi


# perform actions according to the input parameters processed
if [ $TRAVIS_INSTALL ] ; then

  OS="travis"
  conda_install
  conda_test

else 

  if [ $OS_PKGS ] ; then
     install_os_packages
  fi

  if [ $INSTALL_SCRIPTS ] || [ $INSTALL_DEVEL ] ; then
     conda_install
  fi

  if [ $INSTALL_TEST ] ; then
     conda_test
  fi

  if [ $INSTALL_UPDATE ] ; then
     conda_update
  fi

  if [ $UNINSTALL ] ; then
     uninstall
  fi

fi # if-variables

