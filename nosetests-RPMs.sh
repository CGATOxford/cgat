#!/usr/bin/env bash

INIT_DIR=`pwd`

# libpq
wget http://yum.postgresql.org/9.3/redhat/rhel-6-x86_64/pgdg-sl93-9.3-1.noarch.rpm
rpm -i pgdg-sl93-9.3-1.noarch.rpm 
yum install postgresql93-devel

# GCProfile
yum install glibc.i686
yum install compat-libstdc++-33.i686

# create a new folder to store external tools
mkdir -p $HOME/CGAT/external-tools
cd $HOME/CGAT/external-tools

# wigToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
chmod +x wigToBigWig
PATH=$PATH:$HOME/CGAT/external-tools

# BEDtools
curl -L https://github.com/arq5x/bedtools2/releases/download/v2.18.2/bedtools-2.18.2.tar.gz > bedtools-2.18.2.tar.gz
tar xzvf bedtools-2.18.2.tar.gz
cd bedtools-2.18.2
make
PATH=$PATH:$HOME/CGAT/external-tools/bedtools-2.18.2/bin

# bx-python
cd ..
export C_INCLUDE_PATH=$HOME/CGAT/virtualenv-1.10.1/cgat-venv/lib/python2.7/site-packages/numpy/core/include

# GCProfile
wget http://tubic.tju.edu.cn/GC-Profile/download/GCProfile_LINUX.tar
tar xvf GCProfile_LINUX.tar
cp GCProfile_LINUX/GCProfile .
cp GCProfile_LINUX/gnuplot .

# Set up other environment variables
cd $INIT_DIR
export PYTHONPATH=$PYTHONPATH:$INIT_DIR
source $HOME/CGAT/virtualenv-1.10.1/cgat-venv/bin/activate

# setup.py develop
python setup.py develop

# run nosetests
nosetests -v tests/test_scripts.py >& nosetests.out

