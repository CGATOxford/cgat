##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
gpipe/export_code.py - package all code from a pipeline
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script analyses makefile based pipelines and
copies all the python scripts and libraries needed
for running it into a single package.


Usage
-----

Example::

   python gpipe/export_code.py --help

Type::

   python gpipe/export_code.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import os
import shutil
import tempfile
import CGAT.Experiment as E

# Needs refactoring with CGATPipelines
try:
    from CGATPipelines.Local import getMakefiles, getScripts, getModules
except ImportError:
    pass


USAGE = """python %s
""" % sys.argv[0]

PREAMBL = """################################################################################
#   Gene prediction pipeline
#
#   Copyright (C) 2007 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""

README_INSTALLATION = """
CODE INSTALLATION

Unpack code into a directory of your choice:

        tar -xvf %(dir)s.tgz

Check if all python system libraries are present:

        python %(dir)s/check_gpipe/setup.py

Setup the master makefile:

        python %(dir)s/gpipe/setup.py

This will create a Makefile in the current directory, which will then be used
as the build directory. Use the -d option to use a different directory or call
the script from the build directory. The Makefile contains numerous options to
tune the pipeline and to connect with other software (for example, database and
connection details to the postgres database).
"""

README_GENERIC = """
DISCLAIMER:

This pipeline was developed for the resources and setup we have
here. Though the scripts (python and perl) and makefiles should be
portable thanks to the widespread distribution of perl, python and
make, you might have to dig into code if your setup varies from the
one we have.  """

README_GPIPE = """
SETUP

Documentation can be found at

http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/cgat/pipelines/Gpipe.html

"""


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/export_code.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-t", "--tool-set", dest="tool_set",
                      help="set to use.", type="choice",
                      choices=("geneprediction", "codonbias", "geneprediction",
                               "codeml_benchmark"))

    parser.add_option("--test", dest="test", action="store_true",
                      help="testing mode [%default]")

    parser.add_option("--release", dest="version", type="string",
                      help="the version string [%default]")

    parser.set_defaults(
        tool_set="geneprediction",
        libdirs=("/home/andreas/cgat/", ),
        scriptdirs=("/home/andreas/cgat/",),
        name=None,
        version="0.0.1",
        test=False,
    )

    (options, args) = E.Start(parser)

    source_directory = os.path.realpath(os.path.dirname(sys.argv[0])) + "/"

    if options.tool_set == "geneprediction":
        makefiles = set(("Makefile.gpipe",))
        method = "gpipe"
        readme = README_GPIPE
        if options.name is None:
            options.name = "gpipe"
    elif os.path.exists("%s/Makefile.%s" % (source_directory, options.tool_set)):
        makefiles = set(("Makefile.%s" % options.tool_set,))
        method = options.tool_set
        readme = README_GENERIC
        if options.name is None:
            options.name = options.tool_set
    else:
        raise ValueError("unknown toolset %s" % options.tool_set)

    if options.test:
        tmp_dir = "tmp"
        os.mkdir(tmp_dir)
    else:
        tmp_dir = tempfile.mkdtemp()

    nmissed_makefiles, nmissed_scripts, nmissed_modules = 0, 0, 0

    makefiles = getMakefiles(
        makefiles, os.path.join(source_directory, "makefiles"))
    target_directory = "%s/%s-%s/makefiles/" % (tmp_dir,
                                                options.name, options.version)
    os.makedirs(target_directory)

    for makefile in makefiles:

        src = os.path.join(source_directory, makefile)

        if os.path.exists(src):
            shutil.copy(src, target_directory)
            E.info("makefile added: %s" % makefile)
        else:
            E.warn("makefile missed: %s" % makefile)
            nmissed_makefiles += 1

    scripts = getScripts(makefiles, source_directory)
    target_directory = "%s/%s-%s/" % (tmp_dir, options.name, options.version)
    modules = set()
    for language, script in scripts:

        # remove any path name variables
        script = os.path.basename(script)
        script = script[script.find(")") + 1:]

        for d in (source_directory,) + options.scriptdirs:
            src = os.path.join(source_directory, script)

            if os.path.exists(src):
                shutil.copy(src, target_directory)
                if language == "python":
                    modules.add(script)
                E.info("script added: %s from %s" % (script, dir))
                break
        else:
            E.warn("script missed: %s" % (script))
            nmissed_scripts += 1

    modules, system_modules = getModules(modules,
                                         (source_directory,) +
                                         options.scriptdirs,
                                         options.libdirs)

    for lib in modules:

        # remove any path name variables
        for dir in (source_directory,) + options.scriptdirs + options.libdirs:
            if os.path.exists(dir + lib):
                shutil.copy(dir + lib,
                            target_directory + lib)

                if options.loglevel >= 1:
                    options.stdout.write(
                        "module added: %s from %s\n" % (lib, dir))
                break
        else:
            if options.loglevel >= 1:
                options.stdout.write("module missed: %s\n" % (lib))
            nmissed_scripts += 1

    # Writing setup check file
    outfile = open(target_directory + "check_gpipe/setup.py", "w")
    for module in system_modules:
        outfile.write("import %s\n" % module[:-3])

    outfile.write(
        "print 'check successfull - all python modules required are present.'\n")

    outfile.close()

    # Writing setup file.
    outfile = open(target_directory + "gpipe/setup.py", "w")
    infile = open(source_directory + "gpipe/setup.py")
    for line in infile:
        if re.search("method = None", line):
            line = line[line.find(method):] + "method = '%s',\n" % method
        outfile.write(line)

    outfile.close()

    # Writing the README file
    params = {'dir': "%s-%s" % (options.name, options.version),
              }

    outfile = open(target_directory + "README", "w")
    outfile.write(PREAMBL + "\n")
    outfile.write(README_INSTALLATION % params + "\n")
    outfile.write(readme % params + "\n")
    outfile.close()

    # Wrap the whole thing up
    os.system("tar -C %s -czf %s/%s-%s.tgz %s-%s" % (tmp_dir, tmp_dir, options.name, options.version,
                                                     options.name, options.version, ))

    E.info("nmissed: makefiles=%i, scripts=%i, modules=%i" %
           (nmissed_makefiles, nmissed_scripts, nmissed_modules))

    if not options.test:
        shutil.copy("%s/%s-%s.tgz" %
                    (tmp_dir, options.name, options.version), ".")
        shutil.rmtree(tmp_dir)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
