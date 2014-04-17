.. _CGATInstallationDocker:

======================
Docker: CGAT Container
======================

Docker_ is an open-source software to easily create lightweight, 
portable, self-sufficient containers from any application. 

CGAT provides a Docker Container with the Code Collection installed.
To use the CGAT Container you need to install Docker by following
the instructions here:

http://docs.docker.io/en/latest/installation/

Once Docker is available in your system you can get the CGAT Container
by doing::

  docker pull cgat/scripts

When you install the CGAT Code Collection via ``pip``, as explained 
in the :ref:`CGATInstallation`, you use::

  cgat --help

to get CLI help, and::

  cgat <script>

to run ``<script>``. Now, with the CGAT Container you do::

  docker run cgat/scripts

to get CLI help, and::

  docker run cgat/scripts <script>

to run ``<script>``

CGAT Container has been tested with Docker 0.9 and up to this version
the Docker client needs root permissions to run. This is a problem for 
the end users willing to use the CGAT Container in a workstation where
they do not have root privileges. However, it is expected that the Docker
client will eventually find a solution to this problem in future versions.

.. _Docker: https://www.docker.io

