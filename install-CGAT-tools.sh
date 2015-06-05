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

RESULT=0

if [ -z "$CGAT_HOME" ] ; then

   if [ -d "$HOME/cgat-install/conda-install" ] ; then
      UNINSTALL_DIR="$HOME/cgat-install"
      RESULT=1
   fi


else

   if [ -d "$CGAT_HOME/conda-install" ] ; then
      UNINSTALL_DIR="$CGAT_HOME"
      RESULT=1
   fi

fi

return $RESULT

}

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

if [ ! -z "$UNINSTALL_DIR" ] ; then

   echo
   echo " An installation of the CGAT code was found in: $UNINSTALL_DIR"
   echo " Please use --location to install CGAT code in a different location "
   echo " or uninstall the current version before proceeding."
   echo
   echo " Installation is aborted."
   echo
   exit 1

fi

# check installation type
if [ "$OS" != "travis" ] ; then
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
   else
      echo
      echo " Wrong installation type! "
      echo " Installation aborted. "
      echo
      exit 1
   fi # if install type
fi # if travis install

# check installation target
if [ "$OS" == "travis" ] ; then

   CONDA_INSTALL_TYPE="cgat-devel"
   export CONDA_INSTALL_DIR=$TRAVIS_BUILD_DIR/conda-install
   cd $TRAVIS_BUILD_DIR

else

   # Go to CGAT_HOME to continue with installation
   if [ -z "$CGAT_HOME" ] ; then
      # install in default location
      export CGAT_HOME=$HOME/cgat-install
      mkdir -p $CGAT_HOME
   fi

   export CONDA_INSTALL_DIR=$CGAT_HOME/conda-install
   cd $CGAT_HOME

fi # if-OS

# download and install conda
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p $CONDA_INSTALL_DIR
export PATH="$CONDA_INSTALL_DIR/bin:$PATH"
hash -r

# install cgat environment
conda update -q conda --yes
conda info -a
conda create -q -n $CONDA_INSTALL_TYPE $CONDA_INSTALL_TYPE gcc=4.8.3 --override-channels --channel https://conda.binstar.org/cgat --channel defaults --channel https://conda.binstar.org/r --channel https://conda.binstar.org/asmeurer --yes

# if installation is 'devel' (outside of travis), checkout latest version from github
if [ "$OS" != "travis" ] ; then

   if [ "$INSTALL_DEVEL" == "1" ] ; then

      # get latest version from Git Hub
      git clone https://github.com/CGATOxford/cgat.git $CGAT_HOME/cgat-code
      cd $CGAT_HOME/cgat-code

      # activate cgat environment
      source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE

      # Set up other environment variables
      export CFLAGS=$CFLAGS" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
      export CPATH=$CPATH" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
      export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
      export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
      export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib/R/lib

      # Python preparation
      python setup.py develop

      if [ ! $? -eq 0 ] ; then
         echo
         echo " There was a problem doing: 'python setup.py develop' "
         echo " Installation did not finish properly. "
         echo 
         echo " Please submit this issue via Git Hub: "
         echo " https://github.com/CGATOxford/cgat/issue "
         echo 
      fi # if-$?
   fi # if INSTALL_DEVEL

   clear
   echo 
   echo " The CGAT code was successfully installed!"
   echo
   echo " To activate the CGAT environment type: "
   echo " $ source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE"
   [ "$INSTALL_SCRIPTS" == "1" ] && echo " cgat --help"
   echo
   echo " To deactivate the environment, use:"
   echo " $ source deactivate"
   echo

fi

} # conda install


# test code with conda install
conda_test() {

# setup environment and run tests
if [ "$OS" == "travis" ] ; then

   CONDA_INSTALL_DIR=$TRAVIS_BUILD_DIR/conda-install

   # activate cgat environment
   if [ "$INSTALL_LITE" == "1" ] ; then
      source $CONDA_INSTALL_DIR/bin/activate cgat-devel-lite
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/cgat-devel-lite/lib/R/lib
   elif [ "$INSTALL_FULL" == "1" ] ; then
      source $CONDA_INSTALL_DIR/bin/activate cgat-devel
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/cgat-devel/lib/R/lib
   else
      source $CONDA_INSTALL_DIR/bin/activate cgat-devel
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_INSTALL_DIR/envs/cgat-devel/lib/R/lib
   fi

   # configure environment
   export CFLAGS=$CFLAGS" -O0 -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
   export CPATH=$CPATH" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
   export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
   export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
   export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib/R/lib

   cd $TRAVIS_BUILD_DIR
   python setup.py develop

   # run nosetests
   if [ "$TEST_IMPORT" == "1" ] ; then
      nosetests -v tests/test_import.py ;
   elif [ "$TEST_STYLE" == "1" ] ; then
      nosetests -v tests/test_style.py ;
   elif [ "$TEST_CMDLINE" == "1" ] ; then
      echo -e "restrict:\n    manifest:\n" > tests/_test_commandline.yaml
      nosetests -v tests/test_commandline.py ;
   else
      nosetests -v tests/test_scripts.py ;
   fi

else

   detect_cgat_installation

   if [ -z "$UNINSTALL_DIR" ] ; then

      echo
      echo " The location of the CGAT code was not found. "
      echo " Please install before testing."
      echo
      exit 1

   fi
   
   if [ -z "$CGAT_HOME" ] ; then
      # default location for cgat installation
      export CGAT_HOME=$HOME/cgat-install
   fi

   CONDA_INSTALL_DIR=$CGAT_HOME/conda-install
   CONDA_INSTALL_TYPE=`ls $CONDA_INSTALL_DIR/envs`
   cd $CGAT_HOME

   if [ "$CONDA_INSTALL_TYPE" == "cgat-scripts-lite" ] || [ "$CONDA_INSTALL_TYPE" == "cgat-scripts" ] ; then
      wget https://github.com/CGATOxford/cgat/archive/v0.2.3.tar.gz
      tar xzf v0.2.3.tar.gz 
      rm v0.2.3.tar.gz
      cd cgat-0.2.3/
      OUTPUT_DIR=`pwd`
   elif [ "$CONDA_INSTALL_TYPE" == "cgat-devel-lite" ] || [ "$CONDA_INSTALL_TYPE" == "cgat-devel" ] ; then
      cd cgat-code
      OUTPUT_DIR=`pwd`
   else
      echo
      echo " There was an error running the tests. "
      echo " Execution aborted. "
      echo
      exit 1
   fi

   # activate cgat environment
   source $CONDA_INSTALL_DIR/bin/activate $CONDA_INSTALL_TYPE

   # Set up other environment variables
   export CFLAGS=$CFLAGS" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
   export CPATH=$CPATH" -I/usr/include/x86_64-linux-gnu -I$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include -L/usr/lib/x86_64-linux-gnu -L$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib"
   export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
   export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/include
   export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib:$CONDA_INSTALL_DIR/envs/$CONDA_INSTALL_TYPE/lib/R/lib

   # Python preparation
   python setup.py develop

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
   fi

fi # if-OS

} # conda_test


# update conda installation
conda_update() {

detect_cgat_installation

if [ -z "$UNINSTALL_DIR" ] ; then

   echo
   echo " The location of the CGAT code was not found. "
   echo " Please install before updating."
   echo
   exit 1

fi

if [ -z "$CGAT_HOME" ] ; then
   # default location for cgat installation
   export CGAT_HOME=$HOME/cgat-install
fi

CONDA_INSTALL_DIR=$CGAT_HOME/conda-install
CONDA_INSTALL_TYPE=`ls $CONDA_INSTALL_DIR/envs`
cd $CGAT_HOME

source $CGAT_HOME/conda-install/bin/activate $CONDA_INSTALL_TYPE
conda update --all

if [ ! $? -eq 0 ] ; then

   echo
   echo " There was a problem updating the installation. "
   echo 
   echo " Please submit this issue via Git Hub: "
   echo " https://github.com/CGATOxford/cgat/issue "
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
echo " ./install-CGAT-tools.sh --cgat-scripts [--location <path/to/folder>]"
echo
echo " The default location is: $HOME/cgat-install"
echo
echo " Otherwise, if you prefer to use the scripts and the pipelines altogether instead, type:"
echo " ./install-CGAT-tools.sh --cgat-devel [--location <path/to/folder>]"
echo
echo " Both installations create a new Conda environment ready to run the CGAT code."
echo " On the other hand, if you are looking for other advanced installation options please visit:"
echo " https://www.cgat.org/downloads/public/cgat/documentation/CGATInstallation.html"
echo 
echo " To test the installation:"
echo " ./install-CGAT-tools.sh --test [--location </full/path/to/folder/without/trailing/slash>]"
echo
echo " To update the Conda packages:"
echo " ./install-CGAT-tools.sh --update [--location <path/to/folder>]"
echo 
echo " To uninstall the CGAT code:"
echo " ./install-CGAT-tools.sh --uninstall [--location <path/to/folder>]"
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
# variable to store input parameters
INPUT_ARGS=$(getopt -n "$0" -o h0123456789: --long "help,
                                                  travis,
                                                  install-os-packages,
                                                  cgat-scripts,
                                                  cgat-devel,
                                                  lite,
                                                  full,
                                                  test,
                                                  update,
						  uninstall,
                                                  location:"  -- "$@")
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

  else

    help_message

  fi # if-args
  

done # while-loop

# sanity checks
if [ "$INSTALL_LITE" == "1" ] && [ "$INSTALL_FULL" == "1" ] ; then

   echo 
   echo " Incorrect input arguments: mixing --full and --lite options is not permitted."
   echo " Installation aborted. Please run -h option."
   echo
   exit 1

elif [ "$INSTALL_SCRIPTS" == "1" ] && [ "$INSTALL_DEVEL" == "1" ] ; then

   echo
   echo " Incorrect input arguments: mixing --cgat-scripts and --cgat-devel is not permitted."
   echo " Installation aborted. Please run -h option."
   echo
   exit 1

fi


# perform actions according to the input parameters processed
if [ "$TRAVIS_INSTALL" == "1" ] ; then

  OS="travis"
  conda_install
  conda_test

else 

  if [ "$OS_PKGS" == "1" ] ; then
     install_os_packages
  fi

  if [ "$INSTALL_SCRIPTS" == "1" ] || [ "$INSTALL_DEVEL" == "1" ] ; then
     conda_install
  fi

  if [ "$INSTALL_TEST" == "1" ] ; then
     conda_test
  fi

  if [ "$INSTALL_UPDATE" == "1" ] ; then
     conda_update
  fi

  if [ "$UNINSTALL" == "1" ] ; then
     uninstall
  fi

fi # if-variables

