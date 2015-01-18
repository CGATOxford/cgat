.. _CGATInstallationDocker:

======================
Docker: CGAT Container
======================

Docker_ is an open-source software to easily create lightweight, 
portable, self-sufficient containers from any application. 

CGAT provides a Docker Container with the Code Collection pre-installed.
To use the CGAT Container you need to install Docker by following
the instructions here:

http://docs.docker.com/installation/

Once Docker is available on your system you can get the CGAT Container
by doing::

  docker pull cgat/scripts

Now you can type::

  docker run cgat/scripts --help

to get command-line help. Then, any script can be run by typing::

  docker run cgat/scripts <script-name>

CGAT Container has been built with Docker 1.3.2 and up to this version
the Docker client needs root permissions to run. This is a problem for 
the end users willing to use the CGAT Container in a workstation where
they do not have root privileges. However, it is expected that the Docker
client will eventually find a solution to this problem in future versions.

.. _Docker: https://www.docker.com

