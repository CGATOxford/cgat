################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2013 Andreas Heger
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
cgat2rdf.py - create rdf description of CGAT script
====================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script creates an rdf description of a CGAT
script.

Optionally, the script outputs also a galaxy xml
description of the scripts' interface.

Usage
-----

Example::

   python cgat2rdf.py bam2stats.py 

Type::

   python cgat2rdf.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import datetime
import collections

from jinja2 import Template

import CGAT.Experiment as E

# handle to original E.Start function
ORIGINAL_START = None

# Parser object collected from child script
PARSER = None

from rdflib import Graph
from rdflib import Namespace
from rdflib.namespace import RDF, RDFS, DCTERMS
from rdflib import Graph, Namespace, Literal, BNode, URIRef
from rdflib.collection import Collection

#DCTerms = Namespace('http://purl.org/dc/terms/')
FOAF = Namespace('http://xmlns.com/foaf/1.1/')
Component = Namespace('http://www.isi.edu/ikcap/Wingse/componentOntology.owl#')
FO = Namespace('http://www.isi.edu/ikcap/Wingse/fileOntology.owl#')
CLP = Namespace('http://www.humgen.nl/climate/ontologies/clp#')

def _e(string):
    return string.replace(' ', '_')

MAP_FORMATS = {
    'tsv' : 'tabular',
    'table' : 'tabular',
    'stats' : 'tabular',
    'csv' : 'tabular',
    }

MAP_TYPE2FORMAT = {
    'gff' : 'gff',
    'gtf' : 'gtf',
    'bam' : 'bam',
    'sam' : 'sam',
    'bigwig': 'bigWig',
    'bed' : 'bed',
    }


class DummyError( Exception ): pass

class Generator:
    '''inspired by:
    https://github.com/zuotian/CLI-mate/blob/master/climate/utils/legacy_parser.py
    '''

    def __init__(self):
        self.graph = Graph()



    def _addTriple( self, s, p, o):
        if type(o) in [BNode, URIRef]:
            self.graph.add((s, p, o))
        elif type(o) is list:
            o_list = BNode()
            self.graph.add((s, p, o_list))
            os = Collection(self.graph, o_list)
            for item in o:
                os.append(Literal(item))
        elif o != '':
            self.graph.add((s, p, Literal(o)))
            
    def _generateStatements(self, subject_node, properties):
        """
        properties = {"from_loc" : ["data_table", "no"]
        "access_location : ["/path/to/somewhere/", "yes"]}
        """
        for (key, values) in properties.items():
            if values[1] == 'no': # not a volatile property.
                a_node = BNode()
                self._addTriple(a_node, RDF['type'], RDF['Statement'])
                self._addTriple(a_node, CLP['relatedTo'], Literal(key))
                self._addTriple(a_node, RDF['subject'], subject_node)
                self._addTriple(a_node, RDF['predicate'], CLP['hasProperty'])
                self._addTriple(a_node, RDF['object'], values[0])


    def _generateDependencies(self, ap_node, dependencies):
        if dependencies:
            for dep in dependencies:
                d_node = BNode()
                self._addTriple(d_node, RDF.type, CLP['dependency'])
                self._addTriple(d_node, CLP['hasDependingItem'], BNode(_e(dep['depending_parameter'])))
                self._addTriple(d_node, CLP['dependingCondition'], dep['depending_condition'])
                self._addTriple(d_node, CLP['hasDependentItem'], ap_node)
                self._addTriple(d_node, CLP['dependentScope'], dep['dependent_scope'])
                self._addTriple(d_node, CLP['effect'], dep['dependent_effect'])

    def _addDict(self, data ):
        '''convert a dictionary to an RDF graph.'''

        t_node = BNode(_e(data['name']))
        self._addTriple(t_node, RDF.type, CLP['CommandLineProgramComponentType'])
        self._addTriple(t_node, DCTERMS['label'], data['name'])
        self._addTriple(t_node, DCTERMS['title'], data['binary'])
        self._addTriple(t_node, DCTERMS['description'], data['description'])
        self._addTriple(t_node, Component['hasVersion'], data['version'])
        self._addTriple(t_node, DCTERMS['comment'], data['help'])
        self._generateStatements(t_node, data['property_bag'])

        r_node = BNode()
        self._addTriple(t_node, Component['hasExecutionRequirements'], r_node)
        self._addTriple(r_node, RDF.type, Component['ExecutionRequirements'])
        self._addTriple(r_node, Component['requiresOperationSystem'], Component['Linux']) # TODO
        if data['interpreter'] != '(binary)':
            self._addTriple(r_node, Component['requiresSoftware'], Component[data['interpreter']])
        if data['grid_access_type'] != '-':
            self._addTriple(r_node, CLP['gridAccessType'], data['grid_access_type'])
            self._addTriple(r_node, Component['gridID'], data['grid_access_location'])
        for req in data['requirements']:
            req_node = BNode()
            self._addTriple(r_node, CLP['requiresSoftware'], req_node)
            self._addTriple(req_node, RDF.type, CLP['Software'])
            self._addTriple(req_node, DCTERMS['title'], req['req_name'])
            self._addTriple(req_node, CLP['gridID'], req['req_location'])
            self._addTriple(req_node, CLP['softwareType'], req['req_type'])

        argument_list = BNode('argument_list')
        self._addTriple(t_node, Component['hasArguments'], argument_list)
        #self._addTriple(argument_list, RDF.type, Component['argumentAndPrefixList'])
        argument_nodes = Collection(self.graph, argument_list)

        input_list = BNode('input_list')
        self._addTriple(t_node, Component['hasInputs'], input_list)
        #self._addTriple(input_list, RDF.type, Component['FileOrCollectionList'])
        input_nodes = Collection(self.graph, input_list)

        output_list = BNode('output_list')
        self._addTriple(t_node, Component['hasOutputs'], output_list)
        #self._addTriple(output_list, RDF.type, Component['FileOrCollectionList'])
        output_nodes = Collection(self.graph, output_list)

        for p in data['parameters']:
            ap_node = BNode(_e(p['name']))
            argument_nodes.append(ap_node)
            self._addTriple(ap_node, RDF.type, Component['ArgumentAndPrefix'])

            a_node = BNode(_e(p['name']) + '_arg')
            self._addTriple(ap_node, Component['hasArgument'], a_node)

            choices = []
            if 'choices' in p and p['choices']:
                choices = map(lambda x: x.strip(), p['choices'].split(','))

            p_type = p['type']
            if p_type == 'integer':
                self._addTriple(a_node, RDF.type, FO['Int'])
                try:
                    self._addTriple(a_node, FO['hasIntValue'], int(p['value']))
                    choices = map(lambda x: int(x), choices)
                except ValueError:
                    pass # do nothing if value is not an integer
            elif p_type == 'float':
                self._addTriple(a_node, RDF.type, FO['Float'])
                try:
                    self._addTriple(a_node, FO['hasFloatValue'], float(p['value']))
                    choices = map(lambda x: float(x), choices)
                except ValueError:
                    pass # do nothing if value is not a float
            elif p_type in ['string', 'select']:
                self._addTriple(a_node, RDF.type, FO['String'])
                self._addTriple(a_node, FO['hasStringValue'], p['value'])
            elif p_type in ['input', 'stdin']:
                self._addTriple(a_node, RDF.type, FO['File'])
                self._addTriple(a_node, DCTERMS['format'], p['format'])
                self._addTriple(a_node, Component['hasValue'], p['value'])
                input_nodes.append(a_node)
            elif p_type in ['output', 'stdout', 'stderr']:
                self._addTriple(a_node, RDF.type, FO['File'])
                self._addTriple(a_node, DCTERMS['format'], p['format'])
                self._addTriple(a_node, Component['hasValue'], p['value'])
                output_nodes.append(a_node)
            else:
                self._addTriple(a_node, Component['hasValue'], p['value'])

            if choices:
                choices = map(lambda x: Literal(x), choices)
                choice_list = BNode(_e(p['name'] + '_choice_list'))
                choice_nodes = Collection(self.graph, choice_list, choices)
                self._addTriple(a_node, CLP['hasValueChoices'], choice_list)

            self._addTriple(ap_node, DCTERMS['title'], p['name'])
            self._addTriple(ap_node, DCTERMS['description'], p['description'])
            self._addTriple(ap_node, RDFS.label, p['label'])
            self._addTriple(ap_node, Component['hasPrefix'], p['arg'])
            self._addTriple(ap_node, CLP['hasAlternativePrefix'], p['arg_long'])
            self._addTriple(ap_node, CLP['order'], int(p['rank']))
            self._addTriple(ap_node, CLP['display'], p['display'])
            self._addTriple(ap_node, CLP['minOccurrence'], p['min_occurrence'])
            self._addTriple(ap_node, CLP['maxOccurrence'], p['max_occurrence'])

            self._generateStatements(ap_node, p['property_bag'])
            self._generateDependencies(ap_node, p['dependencies'])
        #for

    def _addMetaInfo(self, data):
        
        # meta data node about the Interface Generator itself.
        ig_node = URIRef('http://climate.host.url')
        self._addTriple(ig_node, RDF.type, FOAF['Agent'])
        self._addTriple(ig_node, DCTERMS['title'], data['meta_title'])
        self._addTriple(ig_node, DCTERMS['creator'], data['meta_title'])
        self._addTriple(ig_node, DCTERMS['hasVersion'], data['meta_title'])

        m_node = URIRef('')
        self._addTriple(m_node, RDF.type, FOAF['Document'])
        self._addTriple(m_node, DCTERMS['creator'], ig_node)
        self._addTriple(m_node, DCTERMS['created'], datetime.datetime.utcnow())
        self._addTriple(m_node, RDFS['label'], 'RDF Definition of ' + data['name'])

    def serialize(self, data, format='n3'):
        #TODO: current RDFLib doesn't support base url serialization!
        base_uri = "http://www.humgen.nl/climate/ontologies/clp" + _e(data['name']) + '.rdf#'
        Base = Namespace(base_uri)
        self.graph.bind('base', Base)

        self.graph.bind('dcterms', DCTERMS)
        self.graph.bind('foaf', FOAF)
        self.graph.bind('co', Component)
        self.graph.bind('fo', FO)
        self.graph.bind('clp', CLP)

        self._addDict(data)
        self._addMetaInfo(data)
        return self.graph.serialize(format=format)

def LocalStart( parser, **kwargs ):
    '''stub for E.Start - set return_parser argument to true'''
    global PARSER
    PARSER = ORIGINAL_START( parser,
                             return_parser = True,
                             **kwargs
                             )
    raise DummyError()

def getDescription( scriptname, docstring ):
    '''get script description from docstring.'''

    description = scriptname
    for line in docstring.split("\n"):
        if line.startswith( scriptname ):
            description = line[line.index("-")+1:].strip()
            break
        
    return description

def guessFormats( scriptname, docstring ):
    '''guess the input/output format of a script.'''

    input_format, output_format = "tsv", "tsv"

    if "2" in scriptname:
        input_format, output_format = scriptname.split("2")
        
    # map CGAT format names to GALAXY ones
    input_format = MAP_FORMATS.get( input_format, input_format )
    output_format = MAP_FORMATS.get( output_format, output_format )

    return input_format, output_format

def buildParam( **kwargs ):
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
    param['type' ] = "text"
    param['dependencies'] = {}
    param['property_bag'] = {}
    param['arg_long'] = '--long-argument'

    param.update( kwargs )
    return param

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                             usage = globals()["__doc__"] )

    parser.add_option( "-f", "--format", dest="output_format", type="choice",
                       choices = ("rdf", "galaxy"),
                       help = "output format [%default]. ")
    
    parser.set_defaults( output_format = "rdf" )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 1:
        raise ValueError( "expected a script name" )

    script_name = args[0]
    if not script_name.endswith(".py"):
        raise ValueError( "expected a python script ending in '.py'" )

    # call other script
    dirname = os.path.dirname( script_name )
    basename = os.path.basename( script_name)[:-3]

    sys.path.insert(0, dirname )
    module = __import__(basename)

    # start script in order to build the command line parser
    global ORIGINAL_START
    ORIGINAL_START = E.Start
    E.Start = LocalStart
    E.info( "loaded modules %s" % module)
    try:
        module.main( argv = ["--help"])
    except DummyError:
        pass

    # get script's docstring
    docstring = module.__doc__

    # for k in dir(PARSER):
    #     print k, getattr(PARSER, k)
    # for option in PARSER.option_list:
    #     print option, option.type, option.help, option._short_opts, option._long_opts, option.default

    #@prefix clp: <http://www.humgen.nl/climate/ontologies/clp#> .
    #@prefix co: <http://www.isi.edu/ikcap/Wingse/componentOntology.owl#> .
    #@prefix dcterms: <http://purl.org/dc/terms/> .

    # n = Namespace("http://example.org/people/")
    g = Generator()                            

    data = collections.defaultdict( str )

    data['meta_title'] = 'Interface generator for CGAT scripts'
    data['meta_author'] = 'Andreas Heger'
    data['meta_version'] = 0.1

    data['name'] = basename
    data['interpreter'] = 'python'
    data['property_bag'] = {}
    data['description'] = getDescription( basename, docstring )
    data['help'] = docstring
    data['version'] = "1.0"
    data['owner'] = "CGAT"
    data['email'] = "andreas.heger@gmail.com"
    data['binary'] = script_name
    
    input_format, output_format = guessFormats( basename, docstring )

    stdin = {}
    stdin['name'] = 'input_file'
    stdin['ns_name'] = 'input_file'
    stdin['type'] = 'stdin'
    stdin['label'] = 'input file'
    stdin['description'] = 'input file'
    stdin['choices'] = None
    stdin['format'] = input_format
    stdin['rank'] = 1
    stdin['display'] = 'show'
    stdin['min_occurrence'] = 1
    stdin['max_occurrence'] = 1
    stdin['value'] = ""
    stdin['arg'] = "&lt;"
    stdin['arg_long'] = ""
    stdin['property_bag'] = {}
    stdin['dependencies'] = {}

    stdout = {}
    stdout['name'] = 'tsvfile'
    stdout['ns_name'] = 'tsvfile'
    stdout['type'] = 'stdout'
    stdout['label'] = 'table'
    stdout['description'] = 'bam file'
    stdout['choices'] = None
    stdout['format'] = output_format
    stdout['rank'] = 1
    stdout['display'] = 'show'
    stdout['min_occurrence'] = 1
    stdout['max_occurrence'] = 1
    stdout['value'] = ""
    stdout['arg'] = "&gt;"
    stdout['arg_long'] = ""
    stdout['property_bag'] = {}
    stdout['dependencies'] = {}

    data['parameters'] = [stdin, stdout ]

    defaults = PARSER.get_default_values()

    # flag to indicate wether script needs to go through cgat_wrapper.py
    use_wrapper = False

    for option in PARSER.option_list:
        # ignore options added by optparse
        if option.dest == None: continue
        
        # ignore benchmarking options
        if option.dest.startswith("timeit"):continue

        # ignore options related to forcing output
        if "force" in option.dest: continue

        # ignore special options:
        if option.dest in ("output_filename_pattern", ): continue

        # ignore output options
        if option.dest in ("stdin", "stdout", "stdlog", "stderr", "loglevel" ): continue

        param = buildParam()

        # get command line option call (long/short option)
        try:
            param['arg'] = option._short_opts[0]
        except IndexError: pass

        try:
            param['arg_long'] = option._long_opts[0]
        except IndexError: pass

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
                    use_wrapper = True
                    data['parameters'].append( buildParam( name = 'wrapper_bam_file',
                                                           ns_name = 'wrapper_bam_file',
                                                           arg_long = '--wrapper-bam-file',
                                                           label = option.dest,
                                                           type = 'data',
                                                           format = 'bam',
                                                           help = option.help,
                                                           value = getattr( defaults,  option.dest ) ) )

                    data['parameters'].append( buildParam( name = 'wrapper_bam_index',
                                                           ns_name = 'wrapper_bam_index',
                                                           arg_long = '--wrapper-bai-file',
                                                           type = 'data',
                                                           value = '${wrapper_bam_file.metadata.bam_index}',
                                                           display = 'hidden' ) )

                    # use long argument
                    data['parameters'].append( buildParam( name = 'wrapper_bam_option',
                                                           ns_name = 'wrapper_bam_option',
                                                           arg_long = '--wrapper-bam-option',
                                                           value = param['arg_long'],
                                                           display = 'hidden' ) )
                    
                    continue

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
        param['value'] = getattr( defaults,  option.dest )

            
        param['dependencies'] = {}
        param['property_bag'] = {}

        if option.dest == "genome_file":
            param['property_bag'] = { 'from_loc' : 'path',
                                      'loc_id' : 'sam_fa',
                                      'loc_id_filter' : '1'}

        data['parameters'].append( param )



    if options.output_format == "rdf":
        print g.serialize(data, format='turtle')  


    elif options.output_format == "galaxy":

        if use_wrapper:
        
            # add hidden option for wrapper
            param = buildParam( 
                name = 'wrapper-command',
                ns_name = 'wrapper-command',
                display = 'hidden',
                type = 'text',
                value = data['binary'],
                label = 'wrapper',
                description = 'wrapper',
                arg_long = "--wrapper-command" )
            
            data['parameters'].append( param )
        
        # point to wrapper
        data['binary'] = "cgat_galaxy_wrapper.py"

        displayMap = collections.defaultdict( list )

        for param in data['parameters']:
            displayMap[ param['display'] ].append( param ) 

        displayMap['normal'] = displayMap['show']
        
        target = Template( open('/ifs/devel/andreas/cgat/scripts/cgat2rdf/galaxy.xml').read() )
        print target.render( data = data, 
                             displayMap = displayMap,
                             outputs = [stdout] )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
