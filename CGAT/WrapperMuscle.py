'''
WrapperMuscle.py - 
======================================================

:Tags: Python

Code
----

'''
import os
import tempfile
import subprocess
import shutil

from types import *
from CGAT import Experiment as Experiment
from CGAT import Mali as Mali


class Error(Exception):

    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)


class ParsingError(Error):

    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line):
        self.message = message + " at line " + line


class UsageError(Error):

    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class Muscle:

    mExecutable = "muscle -maxmb 10000"

    def __init__(self):

        self.mTempdir = None

    def __del__(self):

        if self.mTempdir and self.mTempdir != "tmp":
            os.system("rm -rf %s" % self.mTempdir)

    # ------------------------------------------------------------------------
    def Run(self, mali,
            tree=None,
            dump=0,
            test=False,
            options={}):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameInput = "input"
        self.mFilenameOutput = "output"

        if test:
            print("# temporary directory is %s" % self.mTempdir)

        mali.writeToFile(open(self.mTempdir + "/" + self.mFilenameInput, "w"),
                         format="fasta")

        statement = " ".join((self.mExecutable,
                              "-in %s" % self.mFilenameInput,
                              "-out %s" % self.mFilenameOutput))

        s = subprocess.Popen(statement,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             cwd=self.mTempdir,
                             close_fds=True)

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise UsageError("Error in running %s \n%s\n%s\nTemporary directory in %s" % (
                self.mExecutable, err, out, self.mTempdir))

        if dump:
            print("# stdout output of %s:\n%s\n######################################" % (
                self.mExecutable, out))

        result = Mali.Mali()

        result.readFromFile(
            open("%s/%s" % (self.mTempdir, self.mFilenameOutput), "r"), format="fasta")

        if not test:
            shutil.rmtree(self.mTempdir)

        return result

if __name__ == "__main__":

    parser = E.OptionParser(
        version="%prog version: $Id: WrapperMuscle.py 2784 2009-09-10 11:41:14Z andreas $")

    parser.add_option("-t", "--filename-input-tree", dest="filename_input_tree", type="string",
                      help="filename with tree information.")
    parser.add_option("-T", "--filename-output-tree", dest="filename_output_tree", type="string",
                      help="output filename with tree information.")
    parser.add_option("-p", "--program", dest="program", type="string",
                      help="program to use.")
    parser.add_option("-o", "--options", dest="options", type="string",
                      help="input options.")

    parser.set_defaults(
        filename_input_tree=None,
        filename_output_tree=None,
        program=None,
        options="",
    )

    (options, args) = Experiment.Start(parser)

    Experiment.Stop()
