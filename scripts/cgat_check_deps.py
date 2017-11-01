'''
cgat_check_deps.py - check whether the software dependencies are on your PATH
=============================================================================

Purpose
-------

.. The goal of this script is to provide a list of third-party command-line
programs used in a Python script given as input, and check whether
they are on your PATH. This is useful to identify dependencies across all
CGAT pipelines and module files.

This script takes the path to a Python script, which is expected to call
command-line programs like we do in CGAT pipelines, i.e.:

   statement = """cmd-1 | cmd-2 | cmd-3""""
   P.run()

Programs called other way (e.g. using subprocess) will not be picked up
by this script.

Usage
-----

.. python cgat_check_deps --pipeline </path/to/pipeline_name.py> [--print-summary]

Example::

   python cgat_check_deps --pipeline CGATPipelines/pipeline_annotations.py

Type::

   cgat cgat_check_deps --help

for command line help.

Command line options
--------------------

'''

import os
import shutil
import sys
import re
import ast
import argparse
import bashlex


# inspired by
# https://docs.python.org/3/library/ast.html#module-ast
# http://bit.ly/2rDf5xu
# http://bit.ly/2r0Uv9t
# really helpful, used astviewer (installed in a conda-env) to inspect examples
# https://github.com/titusjan/astviewer
def is_cgat_statement(node):
    '''
       Auxiliary function to check for cgat statement:
           statement = "command"
    '''

    result = False
    result = type(node) is ast.Assign and \
        hasattr(node, 'targets') and \
        hasattr(node.targets[0], 'id') and \
        node.targets[0].id == "statement"

    return result


def is_cgat_executable(node):
    '''
       Auxiliary function to check for cgat statement:
           executable = "command"
    '''

    result = False
    result = type(node) is ast.Assign and \
        hasattr(node, 'targets') and \
        hasattr(node.targets[0], 'id') and \
        node.targets[0].id == "executable"

    return result


def is_cgat_executable_name(node):
    '''
       Auxiliary function to check for cgat statement:
           executable_name = "command"
    '''

    result = False
    result = type(node) is ast.Assign and \
        hasattr(node, 'targets') and \
        hasattr(node.targets[0], 'id') and \
        node.targets[0].id == "executable_name"

    return result


def is_cgat_append(node):
    '''
       Auxiliary function to check for cgat statement:
           statment.append("command")
    '''

    result = False
    result = type(node) is ast.Expr and \
        hasattr(node, 'value') and \
        hasattr(node.value, 'func') and \
        hasattr(node.value.func, 'value') and \
        hasattr(node.value.func.value, 'id') and \
        node.value.func.value.id == "statement" and \
        hasattr(node.value.func, 'attr') and \
        node.value.func.attr == "append"

    return result


def get_cmd_string(node):
    '''
       Auxiliary function to get commands in the cgat statement:
           statement = "command"
    '''

    result = ""
    if hasattr(node.value, 's'):
        result = node.value.s
    elif hasattr(node.value, 'left') and hasattr(node.value.left, 's'):
        result = node.value.left.s

    return result


def get_append_string(node):
    '''
       Auxiliary function to get commands in the cgat statement:
           statement.append("command")
    '''

    result = ""
    if hasattr(node, 'value') and \
       hasattr(node.value, 'args') and \
       hasattr(node.value.args[0], 's'):
        result = node.value.args[0].s
    elif hasattr(node, 'value') and \
            hasattr(node.value, 'args') and \
            hasattr(node.value.args[0], 'left') and \
            hasattr(node.value.args[0].left, 's'):
        result = node.value.args[0].left.s

    return result


def cleanup_statement(statement):
    '''
       Auxiliary function to cleanup cgat statements
    '''
    # cleanup whitespaces, tabs, and newlines
    result = " ".join(statement.split())
    # cleanup parameter interpolation
    result = re.sub("\%\(\w+\)\w+", "CGATparameter", result)
    return result


# Thanks to: https://github.com/idank/bashlex
# Workflow:
# statement = ''' <whatever> '''
# statement = " ".join(statement.split())
# statement = re.sub("\%\(", "", statement)
# statement = re.sub("\)s", "", statement)
# parts = bashlex.parse(statement)
# commands = []
# get_cmd_names(parts[0], commands)
def get_cmd_names(tree, commands):
    if hasattr(tree, 'parts') and len(tree.parts) == 0:
        return
    else:
        if hasattr(tree, 'kind'):
            if tree.kind == 'command' and hasattr(tree.parts[0], 'word'):
                commands.append(str(tree.parts[0].word))
            if (tree.kind == 'processsubstitution' or tree.kind == 'commandsubstitution') and \
                    hasattr(tree, 'command') and hasattr(tree.command, 'parts') and \
                    hasattr(tree.command.parts[0], 'word'):
                commands.append(str(tree.command.parts[0].word))
        if hasattr(tree, 'parts'):
            for e in tree.parts:
                get_cmd_names(e, commands)
        if hasattr(tree, 'command'):
            for e in tree.command.parts:
                get_cmd_names(e, commands)


def checkDepedencies(pipeline):

    # check existence of pipeline script
    if not os.access(pipeline, os.R_OK):
        raise IOError("Pipeline %s was not found\n" % pipeline)

    if os.path.isdir(pipeline):
        raise IOError("The given input is a folder, and must be a script\n")

    # parse pipeline script
    with open(pipeline) as f:
        tree = ast.parse(f.read())

    # list to store all statements = ''' <commands> '''
    statements = []

    # inspired by
    # https://docs.python.org/3/library/ast.html#module-ast
    # http://bit.ly/2rDf5xu
    # http://bit.ly/2r0Uv9t
    # really helpful, used astviewer (installed in a conda-env) to inspect examples
    # https://github.com/titusjan/astviewer
    for node in ast.walk(tree):
        statement = ""
        if is_cgat_statement(node) or \
           is_cgat_executable(node) or \
           is_cgat_executable_name(node):

            statement = get_cmd_string(node)

        elif is_cgat_append(node):
            statement = get_append_string(node)

        if len(statement) > 0 and not statement.startswith(' -'):
            #print(statement)
            statement = cleanup_statement(statement)
            statements.append(statement)

    # dictionary where:
    # key = program name
    # value = number of times it has been called
    deps = {}

    # set of names that are not proper deps
    exceptions = ['create',
                  'drop',
                  'select',
                  'attach',
                  'insert',
                  'module',
                  '%s',
                  'tool',
                  'cmd-farm',
                  'cmd-sql',
                  'cmd_extract',
                  'cmds',
                  'compress',
                  'conda_env',
                  'filter_cmd',
                  'load_statement',
                  'match_executable',
                  'rscript',
                  'paired_processing',
                  'executable',
                  'transform',
                  'uncompress',
                  'execglam2',
                  'execglam2scan',
                  'CGATparameter',
                  'checkpoint',
                  'for']

    for statement in statements:
        # use bashlex to parse statements
        commands = []
        try:
            #print(statement)
            parts = bashlex.parse(statement)
            get_cmd_names(parts[0], commands)
        except bashlex.errors.ParsingError:
            pass

        for command in commands:
            #print(command)
            if command.lower() not in exceptions:
                if command not in deps:
                    deps[command] = 1
                else:
                    deps[command] += 1

    # list of unmet dependencies
    check_path_failures = []

    # print dictionary ordered by value
    for k in sorted(deps, key=deps.get, reverse=True):
        if shutil.which(k) is None:
            check_path_failures.append(k)

    return deps, check_path_failures


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if (sys.version_info < (3, 0, 0)):
        raise OSError("This script is Python 3 only")
        sys.exit(-1)

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = argparse.ArgumentParser(description='Get 3rd party dependencies.')

    parser.add_argument("pipeline", help="Path to CGAT pipeline or module")

    parser.add_argument("-s", "--print-summary", dest="summary",
                        action="store_true", default=False,
                        help="Print how many times a program is used")

    options = parser.parse_args()

    # get dependencies dependencies
    deps, check_path_failures = checkDepedencies(options.pipeline)

    # print info about dependencies
    if len(deps) == 0:
        print('\nNo dependencies found.\n')
    else:
        # print dictionary ordered by value
        if options.summary:
            for k in sorted(deps, key=deps.get, reverse=True):
                print('\nProgram: {0!s} used {1} time(s)'.format(k, deps[k]))

        n_failures = len(check_path_failures)
        if n_failures == 0:
            print('\nCongratulations! All required programs are available on your PATH\n')
        else:
            print('\nThe following programs are not on your PATH')
            for p in check_path_failures:
                print('\n{0!s}'.format(p))
            print


if __name__ == "__main__":
    sys.exit(main(sys.argv))
