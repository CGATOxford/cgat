.. _CGATInstallationVirtualBoxUbuntu:

================================
Virtualbox: CGAT Virtual Machine
================================

CGAT also provides a virtual machine with the CGAT Code 
Collection installed in Ubuntu 12.04 LTS. This virtual 
machine has been created with VirtualBox_. If you do not
have VirtualBox installed, please go to the official page
to download and install a copy:

https://www.virtualbox.org/wiki/Downloads

You also need to download the CGAT virtual machine from:

http://www.cgat.org/downloads/cgat-vm.vmdk

To create a new CGAT virtual machine open VirtualBox 
and go through these steps:

1. Click ``New`` and type:

   - Name: cgat-vm

   - Type: Linux

   - Version: Ubuntu (64 bit)

2. RAM memory:

   - Select: ``1024 MB`` of RAM memory at least.

3. Hard drive

   - Select the option: ``Use an existing virtual hard drive file``

   - Click on the dialog and open the file ``cgat-vm.vmdk`` downloaded previously.

The CGAT virtual machine has been created. Power on this 
machine by clicking ``Start``. You may now want to launch
the ``terminal`` and start using the CGAT Code Collection.
To begin, type::

        cgat --help

In case you need root access to this virtual machine, do::

        sudo su

where the password is ``cgat``.

.. _VirtualBox: https://www.virtualbox.org
