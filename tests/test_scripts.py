"""test_scripts.py
==================

nose test script for CGAT scripts

This script builds test cases from directories in the :file:`tests`
subdirectories. Test data for scripts are contained in a directory
with the name of the script and described in a :file:`tests.yaml`
within that directory.

To permit the parallelization of running tests, tests can be run in
chunks (tasks). The script will look for the following environment variables:

CGAT_TASK_ID
   The starting number of this task starting at 1

CGAT_TASK_STEPSIZE
   The number of tests to run within a chunk

"""

from __future__ import print_function

import subprocess
import tempfile
import os
import shutil
import re
import glob
import gzip
import yaml
import time
import hashlib
import sys
import platform

import TestUtils

from nose.tools import ok_

PYTHON_VERSION = platform.python_version()
IS_PY3 = sys.version_info.major >= 3

TRAVIS = os.environ.get("TRAVIS", False) == "true"
JENKINS = os.environ.get("USER", "") == "jenkins"

SUBDIRS = ("gpipe", "optic")

# Setup logging
LOGFILE = open("test_scripts.log", "a")
DEBUG = os.environ.get("CGAT_DEBUG", False)


def check_main(script):
    '''test is if a script can be imported and has a main function.
    '''

    # The following tries importing, but I ran into problems - thus simply
    # do a textual check for now
    # path, basename = os.path.split(script)
    # pyxfile = os.path.join( path, "_") + basename + "x"
    # ignore script with pyximport for now, something does not work
    # if not os.path.exists( pyxfile ):
    #     with warnings.catch_warnings() as w:
    #         warnings.simplefilter("ignore")
    #         (file, pathname, description) =
    #                imp.find_module( basename[:-3], [path,])
    #         module = imp.load_module( basename, file, pathname, description)
    #     ok_( "main" in dir(module), "no main function" )

    # subsitute gpipe and other subdirectories.
    for s in SUBDIRS:
        script = re.sub("%s_" % s, "%s/" % s, script)

    # check for text match
    ok_([x for x in open(script) if x.startswith("def main(")],
        "no main function")


def compute_checksum(filename):
    '''return md5 checksum of file.'''
    return hashlib.md5(open(filename, 'rb').read()).hexdigest()

#########################################
# List of tests to perform.
#########################################
# The fields are:


def check_script(test_name, script, stdin,
                 options, outputs,
                 references,
                 working_dir,
                 current_dir):
    '''check script.

    # 1. Name of the script
    # 2. Filename to use as stdin
    # 3. Option string
    # 4. List of output files to collect
    # 5. List of reference files
    working_dir - directory of test data
    '''
    tmpdir = tempfile.mkdtemp()

    t1 = time.time()

    stdout = os.path.join(tmpdir, 'stdout')
    if stdin:
        if stdin.endswith(".gz"):
            # zcat on osX requires .Z suffix
            stdin = 'gunzip < %s/%s |' % (os.path.abspath(working_dir), stdin)
        else:
            stdin = 'cat %s/%s |' % (os.path.abspath(working_dir), stdin)
    else:
        stdin = ""

    if options:
        options = re.sub("%TMP%", tmpdir, options)
        options = re.sub("<TMP>", tmpdir, options)
        options = re.sub("%DIR%", os.path.abspath(working_dir), options)
        options = re.sub("<DIR>", os.path.abspath(working_dir), options)
    else:
        options = ""

    options = re.sub("\n", "", options)
    short_name = os.path.basename(script)[:-3]

    # use /bin/bash in order to enable "<( )" syntax in shells
    statement = ("/bin/bash -c '%(stdin)s cgat %(short_name)s "
                 " %(options)s"
                 " > %(stdout)s'") % locals()

    process = subprocess.Popen(statement,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=tmpdir)

    if DEBUG:
        print("tmpdir={}".format(tmpdir), end=" ")

    process_stdout, process_stderr = process.communicate()

    fail = False
    msg = ""

    if process.returncode != 0:
        fail = True
        msg = "error in statement: %s; stderr=%s" %\
              (statement, process_stderr)

    # for version tests, do not compare output
    if test_name == "version":
        pass
    elif not fail:
        # compare line by line, ignoring comments
        for output, reference in zip(outputs, references):
            if output == "stdout":
                output = stdout
            elif output.startswith("<DIR>/") or \
                    output.startswith("%DIR%/"):
                output = os.path.join(working_dir, output[6:])
            else:
                output = os.path.join(tmpdir, output)

            if not os.path.exists(output):
                fail = True
                msg = "output file '%s'  does not exist: %s" %\
                      (output, statement)

            reference = os.path.join(working_dir, reference)
            if not fail and not os.path.exists(reference):
                fail = True
                msg = "reference file '%s' does not exist (%s): %s" %\
                      (reference, tmpdir, statement)

            if not fail:
                a = _read(output)
                b = _read(reference)
                if a != b:
                    fail = True
                    msg = ("files %s and %s are not the same\n"
                           "%s\nmd5: output=%i, %s reference=%i, %s") %\
                        (output, reference, statement,
                         len(a),
                         compute_checksum(output),
                         len(b),
                         compute_checksum(reference))

                    diffs = []
                    for aa, bb in zip(a, b):
                        if aa != bb:
                            diffs.append((aa, bb))
                            if len(diffs) > 10:
                                break

                    msg += "first 10 differences: {}".format(
                        "\n--\n".join(
                            ["\n".join(map(str, (x)))
                             for x in diffs]))
                    break

    t2 = time.time()
    LOGFILE.write("%s\t%s\t%f\n" % (script,
                                    test_name,
                                    t2-t1))
    LOGFILE.flush()

    # preserve coverage information, this gets stored it tmpdir, but
    # needs to be moved to the current directory to be merged.
    coverage_files = glob.glob(os.path.join(tmpdir, ".coverage*"))
    for f in coverage_files:
        shutil.move(os.path.abspath(f),
                    os.path.join(current_dir, os.path.basename(f)))

    if not DEBUG:
        shutil.rmtree(tmpdir)
    ok_(not fail, msg)


def test_scripts():
    '''yield list of scripts to test.'''
    # the current directory
    current_dir = os.getcwd()

    # directory location of tests
    testing_dir = TestUtils.get_tests_directory()

    # directory location of scripts
    scripts_dir = os.path.join(os.path.dirname(testing_dir), "CGAT", "scripts")

    # directories with tests (correspond to script names and
    # hence end in .py)
    test_dirs = glob.glob(os.path.join(testing_dir, "*.py"))

    # the config file
    config_file = os.path.join(testing_dir, "_test_scripts.yaml")

    if os.path.exists(config_file):
        config = yaml.load(open(config_file))
        if config is not None:
            if "restrict" in config and config["restrict"]:
                values = config["restrict"]
                if "glob" in values:
                    test_dirs = os.path.join(testing_dir, values["glob"])
                if "manifest" in values:
                    # take scripts defined in the MANIFEST.in file
                    test_dirs = [x for x in open("MANIFEST.in")
                                 if x.startswith("include CGAT/scripts") and
                                 x.endswith(".py\n")]
                    test_dirs = [re.sub("include\s*CGAT/scripts/", "tests/",
                                        x[:-1]) for x in test_dirs]

                if "regex" in values:
                    rx = re.compile(values["regex"])
                    test_dirs = list(filter(rx.search, test_dirs))

    # ignore those which don't exist as tests (files added through MANIFEST.in,
    # such as version.py, __init__.py, ...
    test_dirs = [x for x in test_dirs if os.path.exists(x)]

    # ignore non-directories
    test_dirs = [x for x in test_dirs if os.path.isdir(x)]

    test_dirs.sort()

    # restrict tests run according to chunk parameters
    starting_test_number = os.getenv('CGAT_TASK_ID', None)
    test_increment = os.getenv('CGAT_TASK_STEPSIZE', None)

    try:
        starting_test_number, test_increment = \
            (int(starting_test_number) - 1,
             int(test_increment))
        test_dirs = test_dirs[starting_test_number:
                              starting_test_number + test_increment]
    except TypeError:
        pass

    for test_script in test_dirs:

        script_name = os.path.basename(test_script)

        check_main.description = os.path.join(script_name, "def_main")
        yield (check_main,
               os.path.abspath(os.path.join(scripts_dir, script_name)))

        fn = '%s/tests.yaml' % test_script
        if not os.path.exists(fn):
            continue

        script_tests = yaml.load(open(fn))

        for test, values in sorted(list(script_tests.items())):
            check_script.description = os.path.join(script_name, test)
            if "skip_python" in values:
                versions = [x.strip() for x in
                            str(values["skip_python"]).split(",")]
                versions = [x for x in versions
                            if PYTHON_VERSION.startswith(x)]
                if len(versions) > 0:
                    continue
            if "skip_travis" in values and TRAVIS:
                continue
            if "skip_jenkins" in values and JENKINS:
                continue
 
            # deal with scripts in subdirectories. These are prefixed
            # by a "<subdir>_" for example: optic_compare_projects.py
            # is optic/compare_procjets.py
            if "_" in script_name:
                parts = script_name.split("_")
                if os.path.exists(os.path.join(
                        "scripts", parts[0], "_".join(parts[1:]))):
                    script_name = os.path.join(parts[0], "_".join(parts[1:]))

            yield(check_script,
                  test,
                  os.path.abspath(os.path.join(scripts_dir, script_name)),
                  values.get('stdin', None),
                  values['options'],
                  values['outputs'],
                  values['references'],
                  test_script,
                  current_dir)


def _read(fn):
    if fn.endswith(".gz"):
        with gzip.open(fn) as inf:
            data = inf.read()
    else:
        with open(fn, "rb") as inf:
            data = inf.read()

    if IS_PY3:
        try:
            data = data.decode("ascii")
        except UnicodeDecodeError:
            return data

    data = [x for x in data.splitlines()
            if not x.startswith("#")]

    return data
