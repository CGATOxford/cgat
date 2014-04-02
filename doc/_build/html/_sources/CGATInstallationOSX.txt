.. _CGATInstallationOSX:

=================
OS X installation
=================

This section describes the installation process on OS X. The installation
below has been tested on a MacBook Pro running OS X 10.8.5.

Begin by installing `homebrew <http://brew.sh/>`_ by following these
`instructions <http://hackercodex.com/guide/mac-osx-mountain-lion-10.8-configuration/>`_

Next, install various packages::

   brew install mysql
   brew install gfortran
   brew install freetype
   brew install Mercurial
   brew install python --with-brewed-openssl

Install the R package from `here <http://cran.r-project.org/bin/macosx/>`_ and 
bedtools_::

   wget http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz
   tar -xvzf http://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz
   cd bedtools-2.17.0/
   make
   cp bin/* /usr/local/bin/

Next install and set up a virtual environment::

   pip install virtualenv

CGAT depends on numerous other python packages not all of which
might install cleanly. Here, we give some more detailed instructions.
Generally we recommend when troubleshooting CGAT installation to do so
within a virtual environment. To create a clean environment, type::

    virtualenv --no-site-packages cgat-python
    source cgat-python/bin/activate

Now, download the list of required packages::

    wget https://raw.github.com/CGATOxford/cgat/master/requires.txt

To install the required basic packages::

    pip install -r requires.txt

Also, bx-python needs to be installed. The current version on pypi is
currently out of date, so to install, do::

    pip install https://bitbucket.org/james_taylor/bx-python/get/tip.tar.bz2

If all of that works, installing the CGAT tools should now be
straight-forward::

    pip install cgat


   






