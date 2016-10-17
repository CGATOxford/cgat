"""Requirements.py - Testing for external dependencies
======================================================

This module contains methods for testing if external dependencies of a
CGAT python module or script are satisfied.

The script relies on a section in a module's doc-string to declare
any external dependencies used directly in the module. A section starts
with the token ``Requirements`` followed by a bullet point list of
programs and executables. The executables map to functions in this
module that obtain a version string for that particular tool. The list
ends at the first blank line::

    Requirements:

    * tophat >= 2.0.13
    * bowtie2 >= 2.2.3
    * bwa >= 0.7.8
    * gsnap >= 2014-01-21
    * star >= 2.3.0e (optional)

If a version is ommitted, the module will only check if a tool or package
is installed, but not its version. If the requirement is followed by the tag
``(optional)``, then a missing dependency will not be flagged as an error.

Dependencies are resolved by a :term:`yaml` formatted input file
containing tool definitions. See :file:`external_dependencies.txt` in
the repository for an example.

Code
----

"""

import re
import os
import sys
import collections
import yaml
from distutils.version import LooseVersion
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
from rpy2.robjects import r as R
import rpy2.rinterface

DEFINITIONS_FILE = "external_dependencies.txt"

DATATYPE = collections.namedtuple(
    "DATATYPE",
    "tool installed_version operation required_version optional status")

REGEX_REQUIREMENTS = re.compile("Requirements:\n+(([*]\s[^\n]*\n)*)\n",
                                re.M)


def readDefinitions(filename):
    '''read definitions from a :term:`yaml` file.'''
    with IOTools.openFile(filename) as f:
        config = yaml.load(f)
        if config is None:
            raise IOError("could not read data from '%s'" % filename)
    return config


def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell.

    Execute the string 'cmd' in a shell with os.popen() and return a 2-tuple
    (status, output).  cmd is actually run as '{ cmd ; } 2>&1', so that the
    returned output will contain output or error messages.  A trailing newline
    is stripped from the output.  The exit status for the command can be
    interpreted according to the rules for the C function wait().  Example:

    >>> import subprocess
    >>> subprocess.getstatusoutput('ls /bin/ls')
    (0, '/bin/ls')
    >>> subprocess.getstatusoutput('cat /bin/junk')
    (256, 'cat: /bin/junk: No such file or directory')
    >>> subprocess.getstatusoutput('/bin/junk')
    (256, 'sh: /bin/junk: not found')
    """
    with os.popen('{ ' + cmd + '; } 2>&1', 'r') as pipe:
        try:
            text = pipe.read()
            sts = pipe.close()
        except:
            process = pipe._proc
            process.kill()
            process.wait()
            raise
    if sts is None:
        sts = 0
    if text[-1:] == '\n':
        text = text[:-1]
    return sts, text


class RequirementChecker(object):
    pass


class ToolChecker(RequirementChecker):

    def __init__(self, tool_definition):

        self.tool_definition = tool_definition

    def isInstalled(self):
        path = IOTools.which(self.tool_definition['executable'])
        if path is None:
            return False
        return True

    def getVersion(self):

        if "version" not in self.tool_definition:
            version = ""
        else:
            version = self.tool_definition["version"]

        retcode, version_text = getstatusoutput(
            self.tool_definition['executable'] + " " + version)

        if "regex" not in self.tool_definition:
            return "?"

        installed_version = re.search(self.tool_definition['regex'],
                                      version_text)

        if not installed_version:
            return None

        installed_version = installed_version.groups()[0]
        return installed_version


@E.cachedfunction
def getRPackageList():
    '''return a dictionary of installed R packages
    mapping to their version.'''

    a = R('''installed.packages(
    fields=c("Package", "Version"))[,c("Package", "Version")]
    ''')
    b = R('''installed.packages(
    fields=c("Package", "Version"))[,c("Version")]
    ''')
    return dict(list(zip(a, b)))


class RPackageChecker(RequirementChecker):
    '''checks if an R-package is installed and can
    be loaded::

       library(...)
       packageVersion(...)

    '''

    def __init__(self, tool_definition):
        self.tool_definition = tool_definition
        self.checkinstalled = False

    def isInstalled(self):
        if self.checkinstalled is False:
            try:
                R('''suppressMessages(library(%s))''' %
                  self.tool_definition['rpackage'])
            except rpy2.rinterface.RRuntimeError:
                return False
        self.checkinstalled = True
        return True

    def getVersion(self):

        self.isInstalled()
        r = R('''packageVersion(%s)''' % self.tool_definition['rpackage'])
        return r


class RPackageQuickChecker(RequirementChecker):
    '''check for an R package by simply testing the list of
    installed packages using the following R snippet::

       packinfo <- installed.packages (fields = c ("Package", "Version"))
       packinfo["graphics",c("Package", "Version")]

    This checker will not check if a package can be loaded
    successfully.

    '''

    def __init__(self, tool_definition):
        self.tool_definition = tool_definition

    def isInstalled(self):
        pl = getRPackageList()
        return self.tool_definition['rpackage'] in pl

    def getVersion(self):
        pl = getRPackageList()
        return pl[self.tool_definition['rpackage']]


def checkRequirementsFromBlock(requirements, counter=None, location=""):
    """check if requirements for a script are satisfied

    This method expects a block with the requirements
    definitions.
    """

    definitions = readDefinitions(DEFINITIONS_FILE)

    if not counter:
        counter = E.Counter()

    results = []

    for line in requirements.split("\n"):
        # skip last empty line
        if not line:
            continue

        is_optional = "(optional)" in line
        if is_optional:
            line = re.sub("\(optional\)", "", line)

        if re.search("[=><]", line):
            tool, required = re.search("[*]\s*(\S*)\s(.*)", line).groups()
            required_op, required_version = re.split("\s+", required.strip())
            assert required_op in ("==", ">=", "<=")
        else:
            tool = re.search("[*]\s*(\S*)", line).groups()[0]
            required_version = "?"

        if tool.lower() not in definitions:
            E.warn("%s: tool %s: no definition found" % (location, tool))
            counter.no_definition += 1
            continue

        tool_definition = definitions[tool.lower()]

        if "executable" in tool_definition:
            checker = ToolChecker(tool_definition)
        elif "rpackage" in tool_definition:
            checker = RPackageQuickChecker(tool_definition)
        else:
            E.warn("%s: tool %s: invalid definition" % (location, tool))
            counter.invalid_definition += 1
            continue

        if not checker.isInstalled():
            E.warn("%s: tool %s: not installed" % (location, tool))
            counter.not_installed += 1
            results.append(DATATYPE._make((tool,
                                           "-",
                                           required_op,
                                           required_version,
                                           is_optional,
                                           False)))
            continue

        if required_version != "?":
            installed_version = checker.getVersion()

            if installed_version is None:
                E.warn("%s: tool %s: could not extract version" %
                       (tool, location))
                counter.no_version += 1
                continue

            if installed_version == "?":
                E.warn("%s: tool %s: installed but can not determine version" %
                       (location, tool))
                counter.unknown_version += 1
                results.append(DATATYPE._make((tool,
                                               installed_version,
                                               required_op,
                                               required_version,
                                               is_optional,
                                               True)))
                continue

            # compare installed version with required version
            comp = LooseVersion(installed_version).__cmp__(required_version)

            ok = (required_op == "==" and comp == 0) or \
                 (required_op == ">=" and comp >= 0) or \
                 (required_op == "<=" and comp <= 0) or \
                 (required_op == ">" and comp == 1) or \
                 (required_op == "<" and comp == -1)
        else:
            required_op = "=="
            installed_version = "?"
            ok = True

        if not ok:
            counter.fail += 1
            counter.version_mismatch += 1
            E.warn("%s: tool %s: %s is not %s %s" %
                   (location, tool,
                    installed_version,
                    required_op,
                    required_version))
        else:
            counter.success += 1

        results.append(DATATYPE._make((tool,
                                       installed_version,
                                       required_op,
                                       required_version,
                                       is_optional,
                                       ok)))
    return results


def checkRequirementsFromFile(filename, counter=None):
    '''check requirements from filename.'''
    assert filename.endswith(".py")
    code = open(filename, "r").read()

    requirements = REGEX_REQUIREMENTS.search(code)
    if requirements is None:
        return []

    requirements = requirements.groups()[0]
    return checkRequirementsFromBlock(
        requirements, counter,
        location=os.path.basename(filename))


def checkRequirementsFromModule(module, counter=None):
    '''check requirements for a specific module.'''
    if module.__doc__ is None:
        return []

    requirements = REGEX_REQUIREMENTS.search(module.__doc__)
    if requirements is None:
        return []

    requirements = requirements.groups()[0]
    return checkRequirementsFromBlock(requirements,
                                      counter,
                                      location=module)


def checkRequirementsFromAllModules():

    all_modules = sys.modules
    counter = E.Counter()
    results = []
    for module in list(sys.modules.keys()):
        if all_modules[module] is not None:
            results.extend(checkRequirementsFromModule(
                all_modules[module],
                counter))
    return counter, results
