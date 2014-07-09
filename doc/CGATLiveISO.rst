.. _CGATLiveISO:

=============
CGAT Live ISO
=============

Here is an ISO image of Ubuntu 12.04 LTS with the CGAT
Code Collection v0.2.1 pre-installed:

https://www.cgat.org/downloads/public/cgat/live-iso/CGAT-Code-Collection-Live-0.2.1-x86_64.iso

You can use this ISO at least in three different ways
to run a live image of the software:

1. :ref:`BurnISO`

2. :ref:`BootableUSB`

3. :ref:`ISOinVBox`

The first two options require a computer to be restarted and
booted either from the DVD or the USB drive whereas the
third option require a virtual machine created within
VirtualBox_ to boot the ISO. After you choose one of the options
you will be able to :ref:`RunLive`.

.. _BurnISO:

Burn the ISO into a DVD
=======================

Open your preferred DVD creator, insert a blank DVD into your
DVD drive and burn the ISO_ image.

The next step is to reboot your computer and ensure that you
boot from the DVD drive. Please be careful and select the
option called `Try Ubuntu` when running the Live DVD so you
do not modify any configuration on your computer. Of course,
it is also possible to install this Live DVD, but it will
delete anything you had on your computer, so it is up to you!

.. _BootableUSB:

Create a bootable USB
=====================

It is possible to convert the ISO_ into a bootable USB image.
Please use the `Startup Disk Creator` software available
in Ubuntu by following these steps:

1. Download the ISO_ image in your Ubuntu desktop.

2. Insert a USB stick with at least 2GB of free space.

3. Open the dash and launch `Startup Disk Creator`.

4. Click `Other` to choose the downloaded ISO file.

5. Select the file and click `Open`.

6. Select the USB stick in the bottom box and click `Make Startup Disk`.

When the process completes you'll be ready to restart your
computer and begin using the live image. Please ensure that
you configure your computer to boot from the USB drive.


.. _ISOinVBox:

Run the ISO in VirtualBox
=========================

If you do not have VirtualBox installed, please go to the
official page to download and install a copy for your operating system:

https://www.virtualbox.org/wiki/Downloads

You need a virtual machine already created in VirtualBox.
If you need help to create one go here:

http://www.virtualbox.org/manual/ch01.html#gui-createvm

When you have a virtual machine that you can use, please follow these steps:

1. Open VirtualBox and select your virtual machine.

2. Click `Settings -> Storage -> CD/DVD Drive`

3. Click on the blank DVD icon, select `Choose a virtual
   CD/DVD disk file`, look for the ISO_ image, select it
   and click `Open`.

4. Select the `Live CD/DVD` option in VirtualBox.

5. Click `OK` and `Start` your Virtual Machine.

.. _RunLive:

Run the Live CGAT Code Collection
=================================

Once the Ubuntu image is loaded either ways please follow
these steps to get the CGAT Code Collection v0.2.1 working:

1. Open a terminal console

2. Load the Python environment::

        source /shared/virtualenv-1.11.6/cgat-venv/bin/activate

3. Run the "cgat" command::

        cgat --help


.. _VirtualBox: https://www.virtualbox.org
.. _ISO: https://www.cgat.org/downloads/public/cgat/livecd/CGAT-Code-Collection-Live-0.2.1-x86_64.iso
