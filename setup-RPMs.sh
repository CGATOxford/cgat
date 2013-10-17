#!/usr/bin/env bash

# RPM prerequisites
yum -y install gcc zlib-devel openssl-devel bzip2-devel gcc-c++ freetype-devel libpng-devel blas atlas lapack gcc-gfortran postgresql-devel R-core-devel readline-devel mysql-devel boost-devel sqlite-devel

# additional configuration for scipy
ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
ln -s /usr/lib64/libatlas.so.3 /usr/lib64/libatlas.so
ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so

