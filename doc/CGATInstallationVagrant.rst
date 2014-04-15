.. _CGATInstallationVagrant:

=================
Vagrant: CGAT Box
=================

Vagrant_ allows you to easily deploy virtual machines
(or *boxes*, in Vagrant's jargon) with your development 
environment installed.

First of all, please download and install Vagrant from:

http://www.vagrantup.com/downloads

Once Vagrant is installed on your end, create a dedicated 
folder to store the CGAT box and go there::

  mkdir $HOME/vagrant-cgat-box
  cd $HOME/vagrant-cgat-box

Now you are ready to initialize Vagrant's environment
by doing::

  vagrant init cgat/precise64

Then you can start up the CGAT box this way::

  vagrant up

Finally, you can log into the CGAT box by typing::

  vagrant ssh

The CGAT box is an Ubuntu 12.04 virtual machine with 
the CGAT Code Collection installed via ``pip``. Therefore, 
the ``cgat`` command is available to use the scripts::

  cgat --help all

We have used Vagrant 1.5 to build the CGAT box.

.. _Vagrant: http://www.vagrantup.com/

