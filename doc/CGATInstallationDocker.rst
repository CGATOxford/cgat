.. _CGATInstallationDocker:

=======================
Docker: CGAT Containers
=======================

Docker_ is an open-source software to easily create lightweight, 
portable, self-sufficient containers from any application. 

Installation
------------

CGAT provides Docker Containers with the Code Collection pre-installed.
To use a CGAT Container you need to install Docker by following
the instructions here:

http://docs.docker.com/installation/

Versions explained
------------------

We have now two Docker containers available via DockerHub_:

* ``cgat/scripts``: comes with the CGAT Code Collection pre-installed (released via PyPi_). This one is considered as the production version of the scripts.

* ``cgat/devel``: comes with the latest development version of the CGAT Code (hosted on GitHub_) and the scripts are still under a testing stage.

Below we explain how to use each of them.

CGAT Code Collection
--------------------

Once Docker is available on your system you can get the CGAT Container with the 
CGAT Code Collection pre-installed by doing::

  docker pull cgat/scripts

Now you can type::

  docker run cgat/scripts --help

to get command-line help. Then, any script can be run by typing::

  docker run cgat/scripts <script-name> -h


Latest development version
--------------------------

On the other hand, if you want to use the latest available version of the scripts
you should run::

  docker pull cgat/devel

Now you can type::

  docker run cgat/devel --help

to get command-line help. Then, any script can be run by typing::

  docker run cgat/devel <script-name> -h


Sharing files/folders between your host computer and the CGAT Container
-----------------------------------------------------------------------

There is a data volume within the CGAT Containers in::

  /shared/data

In order to share data between with a CGAT container, please select
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

CGAT Container has been built with Docker 1.5.0 and up to this version
the Docker client needs root permissions to run. This is a problem for 
the end users willing to use the CGAT Container in a workstation where
they do not have root privileges. However, it is expected that the Docker
client will eventually find a solution to this problem in future versions.

.. _Docker: https://www.docker.com
.. _DockerHub: https://registry.hub.docker.com
.. _PyPi: https://pypi.python.org/pypi/CGAT
.. _GitHub: https://github.com/CGATOxford/cgat
