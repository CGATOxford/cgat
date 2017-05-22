'''cgat2dot.py - create a graph between CGAT scripts
====================================================

:Tags: Python

Purpose
-------

This script creates an rdf description of a CGAT script.

Optionally, the script outputs also a galaxy xml description of the
scripts' interface.

Usage
-----

Example::

   python cgat2dot.py scripts/*.py

Type::

   python cgat2dot.py --help

for command line help.

Documentation
-------------

Command line options
--------------------

'''

import os
import sys
import re
import imp

import CGAT.Experiment as E

BASE_URL = "https://www.cgat.org/downloads/public/cgat/documentation/"

ORIGINAL_START = None

PARSER = None


def _e(string):
    return string.replace(' ', '_')

MAP_FORMATS = {
    'tsv': 'table',
    'table': 'table',
    'stats': 'table',
    'csv': 'table',
}

PRINCIPAL_FORMATS = ('bam',
                     'gff',
                     'gtf',
                     'bed',
                     'wiggle',
                     'fasta',
                     'fastq',
                     'fastqs')

BREAK_FORMATS = {'table': 0}
MAP_TYPE2FORMAT = {
    'gff': 'gff,gtf',
    'gtf': 'gff,gtf',
    'bam': 'bam',
    'sam': 'sam',
    'bigwig': 'bigWig',
    'bed': 'bed',
}

NODE_STYLE_DEFAULT = 'color="#A5BB00",style="filled"'
NODE_STYLE_FORMAT = 'color="#7577B8",style="filled"'

EDGE_STYLE_CONVERSION = 'color="#7577B8",penwidth=2'
EDGE_STYLE_DEFAULT = 'color="#A5BB00",penwidth=1'


class DummyError(Exception):
    pass


def LocalStart(parser, *args, **kwargs):
    '''stub for E.Start - set return_parser argument to true'''
    global PARSER
    PARSER = ORIGINAL_START(parser,
                            return_parser=True,
                            **kwargs
                            )
    raise DummyError()


def getDescription(scriptname, docstring):
    '''get script description from docstring.'''

    description = scriptname
    for line in docstring.split("\n"):
        if line.startswith(scriptname):
            description = line[line.index("-") + 1:].strip()
            break

    return description


def guessFormats(scriptname, docstring):
    '''guess the input/output format of a script.'''

    input_format, output_format = "tsv", "tsv"

    if "2" in scriptname:
        input_format, output_format = scriptname.split("2")

    # map CGAT format names to GALAXY ones
    input_format = MAP_FORMATS.get(input_format, input_format)
    output_format = MAP_FORMATS.get(output_format, output_format)

    return input_format, output_format


def buildParam(**kwargs):
    '''return a parameter with default values.

    Specific fields can be set by providing keyword arguments.
    '''

    param = {}

    param['label'] = "label"
    param['description'] = "description"
    param['rank'] = 1
    param['display'] = 'show'
    param['min_occurrence'] = 0
    param['max_occurrence'] = 1

    # get default value
    param['value'] = "value"
    param['type'] = "text"
    param['dependencies'] = {}
    param['property_bag'] = {}
    param['arg_long'] = '--long-argument'

    param.update(kwargs)
    return param


def processScript(script_name, outfile, options):
    '''process one script.'''

    # call other script
    prefix, suffix = os.path.splitext(script_name)

    dirname = os.path.dirname(script_name)
    basename = os.path.basename(script_name)[:-3]

    if options.src_dir:
        dirname = options.src_dir
        script_name = os.path.join(dirname, basename) + ".py"

    if os.path.exists(prefix + ".pyc"):
        os.remove(prefix + ".pyc")

    pyxfile = os.path.join(dirname, "_") + basename + ".pyx"
    if os.path.exists(pyxfile):
        pass

    try:
        module = imp.load_source(basename, script_name)
    except ImportError as msg:
        E.warn('could not import %s - skipped: %s' % (basename, msg))
        return

    E.info("loaded module %s" % module)

    E.Start = LocalStart
    try:
        module.main(argv=["--help"])
    except TypeError as msg:
        E.warn('could not import %s: %s' % (basename, msg))
        return
    except DummyError:
        pass

    # get script's docstring
    docstring = module.__doc__

    input_format, output_format = guessFormats(basename, docstring)

    if output_format in BREAK_FORMATS:
        nodename = '%s%i' % (output_format, BREAK_FORMATS[output_format])
        outfile.write('%s [label="%s"];\n' %
                      (nodename,
                       output_format))
        BREAK_FORMATS[output_format] += 1
        output_format = nodename

    url = BASE_URL + "scripts/%s.html" % basename

    # Note that URL needs to be uppercase!
    if input_format in PRINCIPAL_FORMATS and \
       output_format in PRINCIPAL_FORMATS:
        edge_style = EDGE_STYLE_CONVERSION
    else:
        edge_style = EDGE_STYLE_DEFAULT
    outfile.write('"%s" -> "%s" [label="%s",URL="%s",%s];\n' %
                  (input_format, output_format, basename, url,
                   edge_style))

    return

    # for k in dir(PARSER):
    #     print k, getattr(PARSER, k)
    # for option in PARSER.option_list:
    # print option, option.type, option.help, option._short_opts,
    # option._long_opts, option.default

    # @prefix clp: <http://www.humgen.nl/climate/ontologies/clp#> .
    # @prefix co: <http://www.isi.edu/ikcap/Wingse/componentOntology.owl#> .
    # @prefix dcterms: <http://purl.org/dc/terms/> .

    defaults = PARSER.get_default_values()

    for option in PARSER.option_list:
        # ignore options added by optparse
        if option.dest is None:
            continue

        # ignore benchmarking options
        if option.dest.startswith("timeit"):
            continue

        # ignore options related to forcing output
        if "force" in option.dest:
            continue

        # ignore some special options:
        # if option.dest in ("output_filename_pattern", ):
        #    continue

        # ignore output options
        if option.dest in ("stdin", "stdout", "stdlog", "stderr", "loglevel"):
            continue

        # remove default from help string
        option.help = re.sub("\[[^\]]*%default[^\]]*\]", "", option.help)

        param = buildParam()

        # get command line option call (long/short option)
        try:
            param['arg'] = option._short_opts[0]
        except IndexError:
            pass

        try:
            param['arg_long'] = option._long_opts[0]
        except IndexError:
            pass

        assert 'arg' in param or 'arg_long' in param

        # print "----------------------------------"
        # print [(x,getattr(option,x)) for x in dir( option )]

        param['name'] = option.dest
        param['ns_name'] = option.dest
        if option.type == "int":
            param['type'] = "integer"
        elif option.type == "float":
            param['type'] = "float"
        elif option.type == "string":
            param['type'] = "text"
            if option.metavar:
                mvar = option.metavar.lower()
                if mvar in MAP_TYPE2FORMAT:
                    param['format'] = MAP_TYPE2FORMAT[mvar]
                    param['type'] = "data"
                if mvar == "bam":
                    pass

        elif option.type == "choice":
            param['type'] = "select"
            param['choices'] = option.choices
            if option.action == "append":
                param['multiple'] = True
        elif option.action.startswith("store"):
            param['type'] = "boolean"
        else:
            raise ValueError("unknown type for %s" % str(option))

        param['label'] = option.dest
        param['description'] = option.help
        param['rank'] = 1
        param['display'] = 'show'
        param['min_occurrence'] = 0
        param['max_occurrence'] = 1

        # get default value
        param['value'] = getattr(defaults,  option.dest)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-f", "--format", dest="output_format", type="choice",
                      choices=("rdf", "galaxy"),
                      help="output format [%default]. ")

    parser.add_option("-l", "--list", dest="filename_list", type="string",
                      help="filename with list of files to export "
                      "[%default]. ")

    parser.add_option("-s", "--source-dir", dest="src_dir", type="string",
                      help="directory to look for scripts [%default]. ")

    parser.add_option("-r", "--input-regex", dest="input_regex", type="string",
                      help="regular expression to extract script name "
                      "[%default]. ")

    parser.add_option("-p", "--output-filename-pattern", dest="output_pattern",
                      type="string",
                      help="pattern to build output filename. Should contain "
                      "an '%s' [%default]. ")

    parser.set_defaults(output_format="rdf",
                        src_dir=None,
                        input_regex=None,
                        output_pattern=None,
                        filename_list=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        E.info("reading script names from stdin")
        for line in options.stdin:
            if line.startswith("#"):
                continue
            args.append(line[:-1].split("\t")[0])

    # start script in order to build the command line parser
    global ORIGINAL_START
    ORIGINAL_START = E.Start

    if options.output_pattern and not options.input_regex:
        raise ValueError(
            "please specify --input-regex when using --output-filename-pattern")

    outfile = options.stdout
    outfile.write("""digraph CGAT {
    size="10,20";
    # scale graph so that there are no overlaps
    overlap=scale;
    splines=True;
\n""")

    # set node format for principal genomic formats
    for format in PRINCIPAL_FORMATS:
        outfile.write('"%s" [shape=box,%s];\n' % (format, NODE_STYLE_FORMAT))

    # general node format
    outfile.write('node [%s];\n' % NODE_STYLE_DEFAULT)

    # go through script to provide edges
    for script_name in args:
        if not script_name.endswith(".py"):
            raise ValueError("expected a python script ending in '.py'")

        E.info("input=%s, output=%s" % (script_name, outfile))
        processScript(script_name, outfile, options)

    outfile.write("}\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
