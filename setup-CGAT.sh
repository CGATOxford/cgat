#!/usr/bin/env bash

# Build Python
cd
mkdir CGAT
wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tgz
tar xzvf Python-2.7.5.tgz
rm Python-2.7.5.tgz
cd Python-2.7.5
./configure --prefix=$HOME/CGAT/Python-2.7.5
make
make install
cd
rm -rf Python-2.7.5

# Create virtual environment
cd
cd CGAT
curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
tar xvfz virtualenv-1.10.1.tar.gz
rm virtualenv-1.10.1.tar.gz
cd virtualenv-1.10.1
$HOME/CGAT/Python-2.7.5/bin/python virtualenv.py cgat-venv
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
pip install -r $HOME/requires.txt
pip install CGAT

# Test CGAT Code Collection
cgat --help

# Print help message
echo
echo
echo "To start using the Python virtual environment with the CGAT code collection, type:"
echo "-> source $HOME/CGAT/virtualenv-1.10.1/cgat-venv/bin/activate"
echo "-> cgat --help"
echo
echo "To finish the Python virtual environment, type:"
echo "->deactivate"
echo
echo
