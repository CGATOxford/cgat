################################################################################
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
#################################################################################
'''
export_code.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python export_code.py --help

Type::

   python export_code.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, string, re, optparse, time, os, shutil, tempfile

import Experiment

USAGE="""python %s
""" % sys.argv[0]

PREAMBL="""################################################################################
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

README_INSTALLATION="""
CODE INSTALLATION

Unpack code into a directory of your choice:

        tar -xvf %(dir)s.tgz

Check if all python system libraries are present:

        python %(dir)s/check_setup.py

Setup the master makefile:

        python %(dir)s/setup.py

This will create a Makefile in the current directory, which will then be used
as the build directory. Use the -d option to use a different directory or call
the script from the build directory. The Makefile contains numerous options to
tune the pipeline and to connect with other software (for example, database and
connection details to the postgres database).
"""

README_GENERIC="""
DISCLAIMER:

This pipeline was developed for the resources and setup we have here. Though the scripts
(python and perl) and makefiles should be portable thanks to the widespread distribution of
perl, python and make, you might have to dig into code if your setup varies from the one we
have.
"""

README_GPIPE="""
SETUP

The gene prediction pipeline requires several input files:

1. peptides.fasta: a set of peptide sequences (transcripts). These are the templates
used for gene finding and transcript prediction. This a standard fasta formatted file.
The identifier is the everything on the first line up to the first white space
character.

2. genome.fasta and genome.idx: An indexed genome file. The index and the file can
be created from a set of fasta formatted input files by the command:

        python %(dir)s/IndexedFasta.py genome files

3. peptides2genes: A map of transcripts to genes. This file contains a tab-separated list
that map transcripts to genes in the input peptide set. The format is:

        gene1   transcript1
        gene1   transcript2
        gene2   transcript3
        ...

4. Check if everything is present. Run

        make check-setup

RUNNING

The pipeline uses makefiles to control script logic. Before executing any make commands,
run

        source setup.csh

to update your paths and other environment variables.

Setup the tables in the databases by running the following command:

        make prepare

This needs to be done only once.

Then run:

        make all

to predict transcripts. Output goes to the file "log" in the working directory.

NOTES:

The pipeline assumes that various bioinformatic tools are in the system path. Requirements are:

   * exonerate
   * seg

Furthermore, the pipeline requires:

   * a postgres database server (http://www.postgresql.org)
   * various python modules (numpy, ...), see check_setup.py
   * alignlib (see http://ekhidna.biocenter.helsinki.fi/downloads/alignlib)) with
     the python extension installed.

We use Sun Grid Engine as job queueing system and assume that for all nodes
the code and data can be reached via the same mount point. All jobs that are run on
the cluster are prefixed by the MAKE variable $(CMD_REMOTE_SUBMIT). You can set this
variable to the empty string to run everything locally or on a mosix cluster:

        make all CMD_REMOTE_SUBMIT=

DISCLAIMER:

This pipeline was developed for the resources and setup we have here. Though the scripts
(python and perl) and makefiles should be portable thanks to the widespread distribution of
perl, python and make, you might have to dig into code if your setup varies from the one we
have.
"""

def getMakefiles( makefiles, source_directory = "", ignore_missing = False):
    """get all makefiles that are included in a set of makefiles.

    Keep the order of inclusion.
    """

    read_makefiles = set()
    output_makefiles = []

    def __getMakefiles( makefiles ):

        new_makefiles = []
        for makefile in makefiles:
            if makefile in read_makefiles: continue
            print makefile
            if os.path.exists( makefile ):
                fn = makefile
            elif os.path.exists( os.path.join( source_directory, makefile) ):
                fn = os.path.join( source_directory, makefile  )
            else:
                if ignore_missing:
                    continue
                else:
                    raise IOError, "could not find %s in %s" % (makefile, source_directory)
                
            output_makefiles.append( fn )
            infile = open(fn, "r")
            
            for line in infile:
                if re.match("include\s+(\S+)", line):
                    fn = re.search("include\s+(\S+)", line).groups()[0]
                    # add explicitely given files
                    if os.path.exists( fn ):
                        new_makefile = fn
                    else:
                        new_makefile = os.path.basename( fn )
                        # remove any path name variables
                        new_makefile = new_makefile[new_makefile.find(")")+1:]
                    if new_makefile not in read_makefiles:
                        new_makefiles.append( new_makefile )
            infile.close()

            read_makefiles.add( makefile )
            
        if new_makefiles:
            __getMakefiles( new_makefiles )

    __getMakefiles( makefiles )
    
    return output_makefiles

def getScripts( makefiles, source_directory ):
    """extract all python and perl scripts from a set of makefiles."""
    
    scripts = set()
    
    for makefile in makefiles:
        for line in open( source_directory + makefile,"r"):        
            if re.search( "python", line):
                try:
                    python_scripts = re.search( "python\s+(\S+.py)", line ).groups()
                except AttributeError:
                    continue
                
                for s in python_scripts: scripts.add( ("python", s) )
    
            if re.search( "perl", line):
                try:
                    perl_scripts = re.search( "perl\s+(\S+.pl)", line ).groups()
                except AttributeError:
                    continue
                
                for s in perl_scripts: scripts.add( ("perl", s) )
    
    return scripts

def getModules( modules, scriptdirs, libdirs):
    """extract all imported libraries (modules) from a set of (python) scripts.

    Libraries that are not found are ignored and assumed to be
    installed systemwide.
    """

    read_modules = set()
    system_modules = set()
    
    def __getModules( modules ):
        
        new_modules = set()

        for lib in modules:
            if lib in read_modules or lib in system_modules: continue
            for x in scriptdirs + libdirs:
                if os.path.exists( x + lib ):
                    for line in open( x + lib,"r"):
                        if re.match("import\s+(\S+)", line):
                            ll = re.search("import\s+(\S+)", line).groups()[0]
                            ll = filter( lambda x: x != "", map(lambda x: x.strip(), ll.split(",")))
                            for l in ll:
                                new_modules.add( l + ".py" )
                                
                    read_modules.add( lib )                                
                    break
            else:
                system_modules.add( lib )

        if new_modules:
            __getModules( new_modules )

    __getModules( modules )

    read_modules = read_modules.difference( modules )
    
    return read_modules, system_modules

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: export_code.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-t", "--tool-set", dest="tool_set",
                      help="set to use.", type="choice",
                      choices=("geneprediction", "codonbias", "geneprediction",
                               "codeml_benchmark" ) )

    parser.set_defaults(
        tool_set = "geneprediction",
        libdirs = ("/net/cpp-group/lib/python/", ),
        scriptdirs = ("/net/cpp-group/scripts/tools/", ),
        name = None,
        version = "0.0.1",
        test = False,
        )

    (options, args) = Experiment.Start( parser )

    source_directory = os.path.realpath(os.path.dirname(sys.argv[0])) + "/"
    
    if options.tool_set == "geneprediction":
        makefiles = set( ("Makefile.gpipe",) )
        method = "gpipe"
        readme = README_GPIPE
        if options.name == None: options.name = "gpipe"
    elif os.path.exists( "%s/Makefile.%s" % (source_directory, options.tool_set) ):
        makefiles = set( ("Makefile.%s" % options.tool_set,) )
        method = options.tool_set
        readme = README_GENERIC
        if options.name == None: options.name = options.tool_set
        
    makefiles = getMakefiles( makefiles, source_directory )

    if options.test:
        tmp_dir = "tmp"
        os.mkdir(tmp_dir)
    else:
        tmp_dir = tempfile.mkdtemp()

    target_directory= "%s/%s-%s/" % ( tmp_dir, options.name, options.version )

    os.mkdir( target_directory )

    nmissed_makefiles, nmissed_scripts, nmissed_modules = 0, 0, 0

    for makefile in makefiles:
        
        if os.path.exists( source_directory + makefile ):
            shutil.copy( source_directory + makefile,
                         target_directory + makefile )
            
            if options.loglevel >= 1:
                options.stdout.write("makefile added: %s\n" % makefile )
        else:
            if options.loglevel >= 1:            
                options.stdout.write("makefile missed: %s\n" % makefile )
            nmissed_makefiles += 1
            
    scripts = getScripts( makefiles, source_directory )
    modules = set()
    for language, script in scripts:

        ## remove any path name variables
        script = os.path.basename( script )        
        script = script[script.find(")")+1:]

        for dir in (source_directory,) + options.scriptdirs:
            if os.path.exists( dir + script ):
                shutil.copy( dir + script,
                             target_directory + script )

                if language == "python":
                    modules.add( script )
                if options.loglevel >= 1:
                    options.stdout.write("script added: %s from %s\n" % (script, dir) )
                break
        else:
            if options.loglevel >= 1:
                options.stdout.write("script missed: %s\n" % (script) )
            nmissed_scripts += 1

    modules, system_modules = getModules( modules,
                                 (source_directory,) + options.scriptdirs,
                                 options.libdirs )

    for lib in modules:

        ## remove any path name variables
        for dir in (source_directory,) + options.scriptdirs + options.libdirs:
            if os.path.exists( dir + lib ):
                shutil.copy( dir + lib,
                             target_directory + lib )

                if options.loglevel >= 1:
                    options.stdout.write("module added: %s from %s\n" % (lib, dir) )
                break
        else:
            if options.loglevel >= 1:
                options.stdout.write("module missed: %s\n" % (lib) )
            nmissed_scripts += 1

    ## Writing setup check file
    outfile = open( target_directory + "check_setup.py", "w")
    for module in system_modules:
        outfile.write("import %s\n" % module[:-3] )

    outfile.write("print 'check successfull - all python modules required are present.'\n")
        
    outfile.close()

    ## Writing setup file.
    outfile = open( target_directory + "setup.py", "w")
    infile = open( source_directory + "setup.py" )
    for line in infile:
        if re.search("method = None", line):
            line = line[line.find(method):] + "method = '%s',\n" % method
        outfile.write(line)

    outfile.close()


    ## Writing the README file
    params = { 'dir':"%s-%s" % (options.name, options.version),
               }
    
    outfile = open( target_directory + "README", "w")
    outfile.write( PREAMBL + "\n" )
    outfile.write( README_INSTALLATION % params + "\n")
    outfile.write( readme % params + "\n" )
    outfile.close()
    
    ## Wrap the whole thing up
    os.system( "tar -C %s -czf %s/%s-%s.tgz %s-%s" % (tmp_dir, tmp_dir, options.name, options.version,
                                                      options.name, options.version, ) )
    
    if options.loglevel >= 1:
        options.stdlog.write("# nmissed: makefiles=%i, scripts=%i, modules=%i\n" % (nmissed_makefiles, nmissed_scripts, nmissed_modules ) )

    if not options.test:
        shutil.copy( "%s/%s-%s.tgz" % (tmp_dir, options.name, options.version) , "." )
        shutil.rmtree( tmp_dir )
        
    Experiment.Stop()
