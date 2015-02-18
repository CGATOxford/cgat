.. _CGATInstallationDocker:

======================
Docker: CGAT Container
======================

Docker_ is an open-source software to easily create lightweight, 
portable, self-sufficient containers from any application. 

Installation
------------

CGAT provides a Docker Container with the Code Collection pre-installed.
To use the CGAT Container you need to install Docker by following
the instructions here:

http://docs.docker.com/installation/

Running CGAT Container
----------------------

Once Docker is available on your system you can get the CGAT Container
by doing::

  docker pull cgat/scripts

Now you can type::

  docker run cgat/scripts --help

to get command-line help. Then, any script can be run by typing::

  docker run cgat/scripts <script-name>

Sharing files/folders between your host computer and the CGAT Container
-----------------------------------------------------------------------

There is a data volume within the CGAT Container in::

  /shared/data

In order to share data between with the CGAT container, please select
a local folder on your host::

  /path/to/shared/folder

with the files/folder you want to share, and then run the container
with the ``-v`` option as follows::

  docker run -v /path/to/shared/folder:/shared/data \
  cgat/scripts <script-name> \
  --stdin=/shared/data/input-file \
  --stdout=/shared/data/output-file

This way you can transfer input/output files between the container and
your docker host.

Other details
-------------

CGAT Container has been built with Docker 1.3.2 and up to this version
the Docker client needs root permissions to run. This is a problem for 
the end users willing to use the CGAT Container in a workstation where
they do not have root privileges. However, it is expected that the Docker
client will eventually find a solution to this problem in future versions.

.. _Docker: https://www.docker.com

