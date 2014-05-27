#!/usr/bin/env bash

# message to display when the OS is not correct
sanity_check_os() {
   echo
   echo " Unsupported operating system "
   echo " " $OS
   echo " Installation aborted "
   echo
   exit 1;
} # sanity_check_os


# function to detect the Operating System
detect_os(){

if [ -f /etc/os-release ]; then

   OS=$(cat /etc/os-release | awk '/VERSION_ID/ {sub("="," "); print $2;}' | sed 's/\"//g' | awk '{sub("\\."," "); print $1;}')
   if [ "$OS" != "12" ] ; then

      echo       
      echo " Ubuntu version not supported "
      echo
      echo " Only Ubuntu 12.x has been tested so far "
      echo 
      exit 1;

   fi

   OS="ubuntu"

elif [ -f /etc/system-release ]; then

   OS=$(cat /etc/system-release | awk ' {print $4;}' | awk '{sub("\\."," "); print $1;}')
   if [ "$OS" != "6" ] ; then
      echo
      echo " Scientific Linux version not supported "
      echo
      echo " Only 6.x Scientific Linux has been tested so far "
      echo
      exit 1;
   fi

   OS="sl"

else

   sanity_check_os

fi
} # detect_os


# function to install operating system dependencies
install_os_packages() {

if [ "$OS" == "ubuntu" -o "$OS" == "travis" ] ; then

   echo
   echo " Installing packages for Ubuntu "
   echo

   sudo apt-get install -y gcc g++ zlib1g-dev libssl-dev libbz2-dev libfreetype6-dev libpng12-dev libblas-dev libatlas-dev liblapack-dev gfortran libpq-dev r-base-dev libreadline-dev libmysqlclient-dev libboost-dev libsqlite3-dev mercurial;

elif [ "$OS" == "sl" ] ; then

   echo 
   echo " Installing packages for Scientific Linux "
   echo

   yum -y install gcc zlib-devel openssl-devel bzip2-devel gcc-c++ freetype-devel libpng-devel blas atlas lapack gcc-gfortran postgresql-devel R-core-devel readline-devel mysql-devel boost-devel sqlite-devel mercurial

   # additional configuration for scipy
   ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
   ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so
   ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so;

else

   sanity_check_os $OS

fi # if-OS
} # install_os_packages

# funcion to install Python dependencies
# by default in $HOME/CGAT
# otherwise, in $CGAT_HOME
install_python_deps() {

if [ "$OS" == "ubuntu" -o "$OS" == "sl" ] ; then

   echo
   echo " Installing Python dependencies for $1 "
   echo

   # Go to CGAT_HOME to continue with installation
   if [ -z "$CGAT_HOME" ] ; then
      # install in default location
      CGAT_HOME=$HOME/CGAT-DEPS
   fi

   # Build Python 2.7.5
   mkdir -p $CGAT_HOME
   cd $CGAT_HOME
   mkdir python_build
   cd python_build
   wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
   tar xzvf Python-2.7.5.tgz
   rm Python-2.7.5.tgz
   cd Python-2.7.5
   ./configure --prefix=$CGAT_HOME/Python-2.7.5
   make
   make install
   cd $CGAT_HOME
   rm -rf python_build/

   # Create virtual environment
   wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
   tar xvfz virtualenv-1.10.1.tar.gz
   rm virtualenv-1.10.1.tar.gz
   cd virtualenv-1.10.1
   $CGAT_HOME/Python-2.7.5/bin/python virtualenv.py cgat-venv
   source cgat-venv/bin/activate

   # Install Python prerequisites
   pip install cython
   pip install numpy
   pip install pysam
   pip install https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2
   pip install biopython
   pip install pybedtools
   pip install matplotlib
   pip install scipy
   pip install -r https://raw.github.com/CGATOxford/cgat/master/requires.txt
   pip install --upgrade setuptools
   pip install CGAT

   # Test CGAT Code Collection
   cgat --help

   # Print help message
   echo
   echo
   echo "To start using the Python virtual environment with the CGAT code collection, type:"
   echo "-> source $CGAT_HOME/virtualenv-1.10.1/cgat-venv/bin/activate"
   echo "-> cgat --help"
   echo
   echo "To finish the Python virtual environment, type:"
   echo "->deactivate"
   echo
   echo ;

elif [ "$OS" == "travis" ] ; then
   # Travis-CI provides a virtualenv with Python 2.7
   echo 
   echo " Installing Python dependencies in travis "
   echo

   # Install Python prerequisites
   pip install cython
   pip install numpy
   pip install pysam
   pip install https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2
   pip install biopython
   pip install pybedtools
   pip install matplotlib
   pip install scipy
   pip install -r https://raw.github.com/CGATOxford/cgat/master/requires.txt
   pip install --upgrade setuptools
   #pip install CGAT ;

else

   sanity_check_os $OS

fi # if-OS
} # install_python_deps

install_nosetests_deps() {

if [ "$OS" == "ubuntu" -o "$OS" == "travis" ] ; then

   # GCProfile
   sudo apt-get install -y libc6-i386 libstdc++5:i386

elif [ "$OS" == "sl" ] ; then

   # libpq
   wget http://yum.postgresql.org/9.3/redhat/rhel-6-x86_64/pgdg-sl93-9.3-1.noarch.rpm
   rpm -i pgdg-sl93-9.3-1.noarch.rpm
   yum install -y postgresql93-devel

   # GCProfile
   yum install -y glibc.i686 compat-libstdc++-33.i686

else

   sanity_check_os

fi # if-OS

} # install_nosetests_deps

# common set of tasks to prepare external dependencies
nosetests_external_deps() {
echo
echo " Running nosetests for $1 "
echo

if [ "$OS" == "travis" ] ; then

   # installation directory in travis
   INIT_DIR_DEPS=`pwd`

   mkdir -p $INIT_DIR_DEPS/external-tools
   cd $INIT_DIR_DEPS/external-tools

elif [ "$OS" == "sl" -o "$OS" == "ubuntu" ] ; then

   # Go to CGAT_HOME to continue with installation
   if [ -z "$CGAT_HOME" ] ; then
      # install in default location
      export CGAT_HOME=$HOME/CGAT-DEPS
   fi

   # create a new folder to store external tools
   mkdir -p $CGAT_HOME/external-tools
   cd $CGAT_HOME/external-tools

else

   sanity_check_os

fi # if-OS

# wigToBigWig
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget --no-check-certificate https://www.cgat.org/downloads/public/external-tools/wigToBigWig
chmod +x wigToBigWig
#PATH=$PATH:$CGAT_HOME/external-tools

# bedGraphToBigWig
wget --no-check-certificate https://www.cgat.org/downloads/public/external-tools/bedGraphToBigWig
chmod +x bedGraphToBigWig

# BEDtools
# curl -L https://github.com/arq5x/bedtools2/releases/download/v2.18.2/bedtools-2.18.2.tar.gz > bedtools-2.18.2.tar.gz
wget --no-check-certificate https://www.cgat.org/downloads/public/external-tools/bedtools-2.18.2.tar.gz
tar xzvf bedtools-2.18.2.tar.gz
rm bedtools-2.18.2.tar.gz
cd bedtools-2.18.2
make
#PATH=$PATH:$CGAT_HOME/external-tools/bedtools-2.18.2/bin

# GCProfile
cd ..
# wget http://tubic.tju.edu.cn/GC-Profile/download/GCProfile_LINUX.tar
wget --no-check-certificate https://www.cgat.org/downloads/public/external-tools/GCProfile_LINUX.tar
tar xvf GCProfile_LINUX.tar
rm GCProfile_LINUX.tar
cp GCProfile_LINUX/GCProfile .
cp GCProfile_LINUX/gnuplot .


if [ "$OS" == "travis" ] ; then
   cd $INIT_DIR_DEPS;
fi

} # nosetests_external_deps


# function to run nosetests
run_nosetests() {

if [ "$OS" == "travis" ] ; then

   # installation where CGAT code collection is cloned
   INIT_DIR=`pwd`

   # GCProfile
   #apt-get install -y libc6-i386 libstdc++5:i386

   # prepare external dependencies
   nosetests_external_deps $OS

   # Set up other environment variables
   cd $INIT_DIR
   export PATH=$PATH:$INIT_DIR/external-tools:$INIT_DIR/external-tools/bedtools-2.18.2/bin
   export PYTHONPATH=$PYTHONPATH:$INIT_DIR

   # bx-python
   export C_INCLUDE_PATH=/home/travis/virtualenv/python2.7/local/lib/python2.7/site-packages/numpy/core/include

   # Python preparation
   python setup.py develop
   python scripts/cgat_rebuild_extensions.py

   # run nosetests
   if [ "$TEST_IMPORT" == "1" ] ; then
      nosetests -v tests/test_import.py ;
   else
      nosetests -v tests/test_scripts.py ;
   fi

elif [ "$OS" == "ubuntu" -o "$OS" == "sl" ] ; then

   # prepare external dependencies
   nosetests_external_deps $OS

   # Go to CGAT_GITHUB to continue with installation
   if [ -z "$CGAT_GITHUB" ] ; then
      # install in default location
      export CGAT_GITHUB=$HOME/CGAT-GITHUB
   fi

   if [ -z "$CGAT_HOME" ] ; then
      # default location for cgat installation
      export CGAT_HOME=$HOME/CGAT-DEPS
   fi

   # clone CGAT repository to run nosetests
   git clone https://github.com/CGATOxford/cgat.git $CGAT_GITHUB
   cd $CGAT_GITHUB

   # Set up other environment variables
   export PATH=$PATH:$CGAT_HOME/external-tools:$CGAT_HOME/external-tools/bedtools-2.18.2/bin
   export PYTHONPATH=$PYTHONPATH:$CGAT_GITHUB
   source $CGAT_HOME/virtualenv-1.10.1/cgat-venv/bin/activate
   
   # bx-python
   export C_INCLUDE_PATH=$CGAT_HOME/virtualenv-1.10.1/cgat-venv/lib/python2.7/site-packages/numpy/core/include

   # Python preparation
   python setup.py develop
   python scripts/cgat_rebuild_extensions.py

   # run tests
   /usr/bin/time -o test_import.time -v nosetests -v tests/test_import.py >& test_import.out
   /usr/bin/time -o test_scripts.time -v nosetests -v tests/test_scripts.py >& test_scripts.out ;

else

   sanity_check_os $OS

fi # if-OS

} # run_nosetests


rerun_nosetests() {

# Go to CGAT_GITHUB to continue with installation
if [ -z "$CGAT_GITHUB" ] ; then
   # default location for cgat installation
   export CGAT_GITHUB=$HOME/CGAT-GITHUB
fi

if [ -z "$CGAT_HOME" ] ; then
   # default location for cgat installation
   export CGAT_HOME=$HOME/CGAT-DEPS
fi

cd $CGAT_GITHUB

# set up environment variables
export PATH=$PATH:$CGAT_HOME/external-tools:$CGAT_HOME/external-tools/bedtools-2.18.2/bin
export PYTHONPATH=$PYTHONPATH:$CGAT_GITHUB
source $CGAT_HOME/virtualenv-1.10.1/cgat-venv/bin/activate
export C_INCLUDE_PATH=$CGAT_HOME/virtualenv-1.10.1/cgat-venv/lib/python2.7/site-packages/numpy/core/include

# Python preparation
python setup.py develop
python scripts/cgat_rebuild_extensions.py

# rerun tests
/usr/bin/time -o test_import.time -v nosetests -v tests/test_import.py >& test_import.out
/usr/bin/time -o test_scripts.time -v nosetests -v tests/test_scripts.py >& test_scripts.out ;

} # rerun_nosetests


# function to display help message
help_message() {
echo
echo " Use this script as follows: "
echo
echo " 1) Become root and install the operating system* packages: "
echo " ./install-CGAT-tools.sh --install-os-packages"
echo
echo " 2) Now, as a normal user (non root), install the Python dependencies** in the default folder ($HOME/CGAT-DEPS): "
echo " ./install-CGAT-tools.sh --install-python-deps"
echo
echo " or specify a custom folder with --cgat-deps-dir option, as follows: "
echo " ./install-CGAT-tools.sh --install-python-deps --cgat-deps-dir /path/to/folder"
echo
echo " At this stage the CGAT Code Collection is ready to go and you do not need further steps. Please type the following for more information:"
if [ -z "$CGAT_HOME" ] ; then
   echo " source $HOME/CGAT-DEPS/virtualenv-1.10.1/cgat-venv/bin/activate"
else
   echo " source $CGAT_HOME/virtualenv-1.10.1/cgat-venv/bin/activate"
fi 
echo " cgat --help "
echo
echo " The CGAT Code Collection tests the software with nosetests. If you are interested in running those, please continue with the following steps:"
echo
echo " 3) Become root to install external tools and set up the environment: "
echo " ./install-CGAT-tools.sh --install-nosetests-deps"
echo
echo " 4) Then, back again as a normal user (non root), run nosetests as follows:"
echo " ./install-CGAT-tools.sh --run-nosetests"
echo 
echo " This will clone the CGAT repository from GitHub to: $HOME/CGAT-GITHUB by default. If you want to change that use --git-hub-dir as follows:"
echo " ./install-CGAT-tools.sh --run-nosetests --git-hub-dir /path/to/folder"
echo
echo " If you wanted to re-run nosetests you can do so by typing:"
echo " ./install-CGAT-tools.sh --rerun-nosetests"
echo
echo " If you installed the CGAT Code Collection with the '--cgat-deps-dir option' you need to include it again to re-run the tests:"
echo " ./install-CGAT-tools.sh --rerun-nosetests --git-hub-dir /path/to/folder"
echo 
echo " NOTES: "
echo " * Supported operating systems: Ubuntu 12.x and Scientific Linux 6.x "
echo " ** An isolated virtual environment will be created to install Python dependencies "
echo
exit 1;
} # help_message

# the script starts here

if [ $# -eq 0 ] ; then

   help_message

fi

# these variables will store the information about input parameters
# travis execution
TRAVIS=
# install operating system's dependencies
OS_PKGS=
# install Python dependencies
PY_PKGS=
# install dependencies to run nosetests
NT_PKGS=
# run nosetests
NT_RUN=
# rerun nosetests
NT_RERUN=
# variable to actually store the input parameters
INPUT_ARGS=$(getopt -n "$0" -o ht1234g:c: --long "help,
                                                  travis,
                                                  install-os-packages,
                                                  install-python-deps,
                                                  install-nosetests-deps,
                                                  run-nosetests,
                                                  rerun-nosetests,
                                                  git-hub-dir:,
                                                  cgat-deps-dir:"  -- "$@")
eval set -- "$INPUT_ARGS"

# process all the input parameters first
while [ "$1" != "--" ]
do

  if [ "$1" == "-h" -o "$1" == "--help" ] ; then

    help_message

  elif [ "$1" == "-t" -o "$1" == "--travis" ] ; then
      
      TRAVIS=1
      shift ;

  elif [ "$1" == "-1" -o "$1" == "--install-os-packages" ] ; then
      
      OS_PKGS=1
      shift ;

  elif [ "$1" == "-2" -o "$1" == "--install-python-deps" ] ; then
      
      PY_PKGS=1
      shift ;

  elif [ "$1" == "-3" -o "$1" == "--install-nosetests-deps" ] ; then

      NT_PKGS=1
      shift ;

  elif [ "$1" == "-4" -o "$1" == "--run-nosetests" ] ; then

      NT_RUN=1
      shift ;

  elif [ "$1" == "-5" -o "$1" == "--rerun-nosetests" ] ; then

      NT_RERUN=1
      shift ;

  elif [ "$1" == "-g" -o "$1" == "--git-hub-dir" ] ; then

      CGAT_GITHUB="$2"
      shift 2 ;

  elif [ "$1" == "-c" -o "$1" == "--cgat-deps-dir" ] ; then

      CGAT_HOME="$2"
      shift 2 ;

  else

    help_message

  fi # if-args
  

done # while-loop

# perform actions according to the input parameters processed
if [ "$TRAVIS" == "1" ] ; then

  OS="travis"
  install_os_packages
  install_python_deps
  install_nosetests_deps
  run_nosetests

else 

  detect_os
  
  if [ "$OS_PKGS" == "1" ] ; then
    install_os_packages
  fi

  if [ "$PY_PKGS" == "1" ] ; then
    install_python_deps
  fi

  if [ "$NT_PKGS" == "1" ] ; then
    install_nosetests_deps
  fi

  if [ "$NT_RUN" == "1" ] ; then
    run_nosetests
  fi

  if [ "$NT_RERUN" == "1" ] ; then
    rerun_nosetests
  fi

fi # if-variables

