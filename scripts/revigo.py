################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
revigo.py - apply semantic clustering to GO output
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

apply semantic clustering ala revigo (:pmid:`21789182`, http://revigo.irb.hr)
to a list of significantly enriched or depleted GO terms.

This script can be used to reduce the redundancy in large lists of
GO results. The principal input is list of GO-terms and associated
P-Values.

Usage
-----

Example::

   python revigo.py

Type::

   python revigo.py --help

for command line help.

Documentation
-------------

In addition to a list of GO terms and their associated p-values (--filename-pvalues)

Additional input files required are:

* a file containing a gene ontology in obo-xml format (--filename-ontology)
* a file mapping gene identifiers to GO terms (--filename-go)

The script will compute a matrix of semantic similarities between GO terms
based on the ``simrel`` measure described in the manuscript by Schlicker et al. (2006)
(pmid:`16776819`).

The script will then remove redundant terms from this matrix following the produce
described in the revigo manuscript (:pmid:`21789182`).

The script will output a file called `cluster.tsv` which contains the clustering results.

It will also output an image in :term:`svg` format of the clustered terms and their
similarity in semantic space. Nodes are coloured by fold change and size is given
by the P-Value.

Additional output files
+++++++++++++++++++++++

Computing the similarity matrix can take same time. The script outputs some intermediate
files that will be re-used on subsequent calls in order to speed up processing. The intermediate 
files are:

    * ancestors.tsv: tab-separated file with all the ancestors of a term
    * termfrequencies.tsv: tab-separated file with all the term frequencies and genes with this term
    * simrel.matrix: the full simrel matrix. This matrix is symmetric.
    * cluster.matrix: the clustered simrel matrix. Only the cluster representatives have been kept.

Quality control
+++++++++++++++

The output of this program is similar but not identical to revigo. The algorithm is not fully described
and the input data are not identical. It seems that simrel values in this script are lower and thus the
clustering is more conservative compared to revigo.

Code
----

This module depends on unreleased code from http://code.google.com/p/gographer for parsing of the
GO graph. The code has been modified to allow parsing of relationship information and record synonyms
and obsolete identifiers. The code has been included into this script until gographer will be released. 

'''

import itertools
import math
import sys
import optparse
import collections
import os
import re
import string

import CGAT.GO as GO
import CGAT.IOTools as IOTools

import CGAT.Experiment as E
import numpy
import networkx
import networkx.algorithms.traversal.depth_first_search as traversal
from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from networkx import topological_sort
import cPickle

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import pylab

class Tokenizer:
    '''from gographer utils.py'''
    def __init__(self):
        self.pattern_punctuation = re.compile('(\!|\.|\+|\?|\;|\,|\:|\&|\_|\~|\:|\-|\;|\$|\@|\#|\%|\*|\^|\`|\|)')
        self.pattern_numbers = re.compile('[0-9]')
        self.pattern_plurals = re.compile('\'s$')
        self.pattern_possessives = re.compile('s\'$')
        self.pattern_blanks = re.compile('[\t, ,\n,\r,\f,\v,]+')
        self.pattern_others = re.compile('(\'|\"|\/|\\\|\<|\>|\{|\}|\[|\]|\(|\))')
        
        
        self.pattern_words_and_blanks = re.compile('[\W\t\n\r\f\v| ]|\|')
        #self.pattern.word = re.compile('[a-z]+')[^a-zA-Z]

        self.leftchar = {}
        self.leftchar[' ']=0
        for i in string.ascii_letters:
            self.leftchar[i]=0
        #for i in range(10):
            #self.leftchar[str(i)]=0
        #print len(self.leftchar)

            
    def keepNumChr(self, words):
        newWords = ''
        for c in words:
            if self.leftchar.has_key(c):
                newWords +=c
        return newWords
    
                                
    def tokenize_word(self,words):
        before = words
        words = words.lower()
        words = self.pattern_plurals.sub('', words)
        words = self.pattern_possessives.sub('', words)
        words = self.keepNumChr(words)
        words = self.pattern_blanks.sub(' ',words)
        
        #words = string.strip(words)    
        return words

        
    def test(self,words):
        for c in words:    
            if not self.leftchar.has_key(c):
                self.leftchar[c] = 1

                
    def returnUnexpChar(self):
        reV = []
        for k in self.leftchar.keys():
            if self.leftchar[k] != 0:
                reV +=[k]
                
        return reV

class PorterStemmer:
    '''from gographer utils.py'''
    def __init__(self):
        """The main part of the stemming algorithm starts here.
        b is a buffer holding a word to be stemmed. The letters are in b[k0],
        b[k0+1] ... ending at b[k]. In fact k0 = 0 in this demo program. k is
        readjusted downwards as the stemming progresses. Zero termination is
        not in fact used in the algorithm.

        Note that only lower case sequences are stemmed. Forcing to lower case
        should be done before stem(...) is called.
        """

        self.b = ""  # buffer for word to be stemmed
        self.k = 0
        self.k0 = 0
        self.j = 0   # j is a general offset into the string

    def cons(self, i):
        """cons(i) is TRUE <=> b[i] is a consonant."""
        if self.b[i] == 'a' or self.b[i] == 'e' or self.b[i] == 'i' or self.b[i] == 'o' or self.b[i] == 'u':
            return 0
        if self.b[i] == 'y':
            if i == self.k0:
                return 1
            else:
                return (not self.cons(i - 1))
        return 1

    def m(self):
        """m() measures the number of consonant sequences between k0 and j.
        if c is a consonant sequence and v a vowel sequence, and <..>
        indicates arbitrary presence,

           <c><v>       gives 0
           <c>vc<v>     gives 1
           <c>vcvc<v>   gives 2
           <c>vcvcvc<v> gives 3
           ....
        """
        n = 0
        i = self.k0
        while 1:
            if i > self.j:
                return n
            if not self.cons(i):
                break
            i = i + 1
        i = i + 1
        while 1:
            while 1:
                if i > self.j:
                    return n
                if self.cons(i):
                    break
                i = i + 1
            i = i + 1
            n = n + 1
            while 1:
                if i > self.j:
                    return n
                if not self.cons(i):
                    break
                i = i + 1
            i = i + 1

    def vowelinstem(self):
        """vowelinstem() is TRUE <=> k0,...j contains a vowel"""
        for i in range(self.k0, self.j + 1):
            if not self.cons(i):
                return 1
        return 0

    def doublec(self, j):
        """doublec(j) is TRUE <=> j,(j-1) contain a double consonant."""
        if j < (self.k0 + 1):
            return 0
        if (self.b[j] != self.b[j-1]):
            return 0
        return self.cons(j)

    def cvc(self, i):
        """cvc(i) is TRUE <=> i-2,i-1,i has the form consonant - vowel - consonant
        and also if the second c is not w,x or y. this is used when trying to
        restore an e at the end of a short  e.g.

           cav(e), lov(e), hop(e), crim(e), but
           snow, box, tray.
        """
        if i < (self.k0 + 2) or not self.cons(i) or self.cons(i-1) or not self.cons(i-2):
            return 0
        ch = self.b[i]
        if ch == 'w' or ch == 'x' or ch == 'y':
            return 0
        return 1

    def ends(self, s):
        """ends(s) is TRUE <=> k0,...k ends with the string s."""
        length = len(s)
        if s[length - 1] != self.b[self.k]: # tiny speed-up
            return 0
        if length > (self.k - self.k0 + 1):
            return 0
        if self.b[self.k-length+1:self.k+1] != s:
            return 0
        self.j = self.k - length
        return 1

    def setto(self, s):
        """setto(s) sets (j+1),...k to the characters in the string s, readjusting k."""
        length = len(s)
        self.b = self.b[:self.j+1] + s + self.b[self.j+length+1:]
        self.k = self.j + length

    def r(self, s):
        """r(s) is used further down."""
        if self.m() > 0:
            self.setto(s)

    def step1ab(self):
        """step1ab() gets rid of plurals and -ed or -ing. e.g.

           caresses  ->  caress
           ponies    ->  poni
           ties      ->  ti
           caress    ->  caress
           cats      ->  cat

           feed      ->  feed
           agreed    ->  agree
           disabled  ->  disable

           matting   ->  mat
           mating    ->  mate
           meeting   ->  meet
           milling   ->  mill
           messing   ->  mess

           meetings  ->  meet
        """
        if self.b[self.k] == 's':
            if self.ends("sses"):
                self.k = self.k - 2
            elif self.ends("ies"):
                self.setto("i")
            elif self.b[self.k - 1] != 's':
                self.k = self.k - 1
        if self.ends("eed"):
            if self.m() > 0:
                self.k = self.k - 1
        elif (self.ends("ed") or self.ends("ing")) and self.vowelinstem():
            self.k = self.j
            if self.ends("at"):   self.setto("ate")
            elif self.ends("bl"): self.setto("ble")
            elif self.ends("iz"): self.setto("ize")
            elif self.doublec(self.k):
                self.k = self.k - 1
                ch = self.b[self.k]
                if ch == 'l' or ch == 's' or ch == 'z':
                    self.k = self.k + 1
            elif (self.m() == 1 and self.cvc(self.k)):
                self.setto("e")

    def step1c(self):
        """step1c() turns terminal y to i when there is another vowel in the stem."""
        if (self.ends("y") and self.vowelinstem()):
            self.b = self.b[:self.k] + 'i' + self.b[self.k+1:]

    def step2(self):
        """step2() maps double suffices to single ones.
        so -ization ( = -ize plus -ation) maps to -ize etc. note that the
        string before the suffix must give m() > 0.
        """
        if self.b[self.k - 1] == 'a':
            if self.ends("ational"):   self.r("ate")
            elif self.ends("tional"):  self.r("tion")
        elif self.b[self.k - 1] == 'c':
            if self.ends("enci"):      self.r("ence")
            elif self.ends("anci"):    self.r("ance")
        elif self.b[self.k - 1] == 'e':
            if self.ends("izer"):      self.r("ize")
        elif self.b[self.k - 1] == 'l':
            if self.ends("bli"):       self.r("ble") # --DEPARTURE--
            # To match the published algorithm, replace this phrase with
            #   if self.ends("abli"):      self.r("able")
            elif self.ends("alli"):    self.r("al")
            elif self.ends("entli"):   self.r("ent")
            elif self.ends("eli"):     self.r("e")
            elif self.ends("ousli"):   self.r("ous")
        elif self.b[self.k - 1] == 'o':
            if self.ends("ization"):   self.r("ize")
            elif self.ends("ation"):   self.r("ate")
            elif self.ends("ator"):    self.r("ate")
        elif self.b[self.k - 1] == 's':
            if self.ends("alism"):     self.r("al")
            elif self.ends("iveness"): self.r("ive")
            elif self.ends("fulness"): self.r("ful")
            elif self.ends("ousness"): self.r("ous")
        elif self.b[self.k - 1] == 't':
            if self.ends("aliti"):     self.r("al")
            elif self.ends("iviti"):   self.r("ive")
            elif self.ends("biliti"):  self.r("ble")
        elif self.b[self.k - 1] == 'g': # --DEPARTURE--
            if self.ends("logi"):      self.r("log")
        # To match the published algorithm, delete this phrase

    def step3(self):
        """step3() dels with -ic-, -full, -ness etc. similar strategy to step2."""
        if self.b[self.k] == 'e':
            if self.ends("icate"):     self.r("ic")
            elif self.ends("ative"):   self.r("")
            elif self.ends("alize"):   self.r("al")
        elif self.b[self.k] == 'i':
            if self.ends("iciti"):     self.r("ic")
        elif self.b[self.k] == 'l':
            if self.ends("ical"):      self.r("ic")
            elif self.ends("ful"):     self.r("")
        elif self.b[self.k] == 's':
            if self.ends("ness"):      self.r("")

    def step4(self):
        """step4() takes off -ant, -ence etc., in context <c>vcvc<v>."""
        if self.b[self.k - 1] == 'a':
            if self.ends("al"): pass
            else: return
        elif self.b[self.k - 1] == 'c':
            if self.ends("ance"): pass
            elif self.ends("ence"): pass
            else: return
        elif self.b[self.k - 1] == 'e':
            if self.ends("er"): pass
            else: return
        elif self.b[self.k - 1] == 'i':
            if self.ends("ic"): pass
            else: return
        elif self.b[self.k - 1] == 'l':
            if self.ends("able"): pass
            elif self.ends("ible"): pass
            else: return
        elif self.b[self.k - 1] == 'n':
            if self.ends("ant"): pass
            elif self.ends("ement"): pass
            elif self.ends("ment"): pass
            elif self.ends("ent"): pass
            else: return
        elif self.b[self.k - 1] == 'o':
            if self.ends("ion") and (self.b[self.j] == 's' or self.b[self.j] == 't'): pass
            elif self.ends("ou"): pass
            # takes care of -ous
            else: return
        elif self.b[self.k - 1] == 's':
            if self.ends("ism"): pass
            else: return
        elif self.b[self.k - 1] == 't':
            if self.ends("ate"): pass
            elif self.ends("iti"): pass
            else: return
        elif self.b[self.k - 1] == 'u':
            if self.ends("ous"): pass
            else: return
        elif self.b[self.k - 1] == 'v':
            if self.ends("ive"): pass
            else: return
        elif self.b[self.k - 1] == 'z':
            if self.ends("ize"): pass
            else: return
        else:
            return
        if self.m() > 1:
            self.k = self.j

    def step5(self):
        """step5() removes a final -e if m() > 1, and changes -ll to -l if
        m() > 1.
        """
        self.j = self.k
        if self.b[self.k] == 'e':
            a = self.m()
            if a > 1 or (a == 1 and not self.cvc(self.k-1)):
                self.k = self.k - 1
        if self.b[self.k] == 'l' and self.doublec(self.k) and self.m() > 1:
            self.k = self.k -1

    def stem(self, p, i, j):
        """In stem(p,i,j), p is a char pointer, and the string to be stemmed
        is from p[i] to p[j] inclusive. Typically i is zero and j is the
        offset to the last character of a string, (p[j+1] == '\0'). The
        stemmer adjusts the characters p[i] ... p[j] and returns the new
        end-point of the string, k. Stemming never increases word length, so
        i <= k <= j. To turn the stemmer into a module, declare 'stem' as
        extern, and delete the remainder of this file.
        """
        # copy the parameters into statics
        self.b = p
        self.k = j
        self.k0 = i
        if self.k <= self.k0 + 1:
            return self.b # --DEPARTURE--

        # With this line, strings of length 1 or 2 don't go through the
        # stemming process, although no mention is made of this in the
        # published algorithm. Remove the line to match the published
        # algorithm.

        self.step1ab()
        self.step1c()
        self.step2()
        self.step3()
        self.step4()
        self.step5()
        return self.b[self.k0:self.k+1]

class GOOboXmlHandler (ContentHandler):
    """ This class handles parsing the OBO XML file from the Gene Ontology Consortium and populate
        a GOGraph object passed to the constructor of this class.  The handler only use the nodes from
        the namespace (a branch of the GO) specificed in the GOGraph object.  Currently, only the edges
        that reflect IS_A relationship between a pair of GO terms are parsed and added into the GOGraph

        Modification by AH: 
        1. also add PART_OF relationships into the graph
        2. return a list of synonyms

        from gographer.GOOboXmlHandler.py
    """
    ## Constructor.
    # @param    goGraph A reference to a GOGraph object to which the parsed GO nodes and edges
    #                   will be added.
    def __init__(self, goGraph):
        self.namespace = goGraph.getNameSpace()
        self.graph = goGraph
        self.inTypedef = False
        self.goNode = GONode()
        self.obsolete = False
        self.relationship_type = None
        self.synonyms = {}
        self.obsoletes = set()

    def startElement(self, name, attributes):
        self.cdata = ""

        if name == "term":
            self.goid = None
            self.name = None
            self.description = None
            self.obsolete = False            
            self.elementNamespace = None
            self.parents = []
            self.goNode = GONode()
            
        elif name == "typedef":
            self.inTypedef = True

    def endElement(self, name):

        # print name, self.obsolete, self.namespace #, self.elementNamespace

        if name == "id":
            self.goid = self.cdata.strip()
            self.goNode.setGOID(self.goid)

        elif name == "namespace":
            self.elementNamespace = self.cdata.strip()
            self.goNode.setNamespace(self.elementNamespace)

        elif name == "is_a" and not self.inTypedef:
            self.parents.append(self.cdata.strip())
            self.goNode.setParents(self.parents)

        elif name == "type":
            self.relationship_type = self.cdata.strip()

        elif name == "relationship" and not self.inTypedef:
            # add relationships into graph - this includes 'part_of' but excludes
            # 'negatively_regulates', 'positively_regulates' and 'regulates' and others
            if self.relationship_type == "part_of":
                self.parents.append(self.cdata.strip())
                self.goNode.setParents(self.parents)
            
        elif name == "is_obsolete" and self.cdata.strip() == "1":
            self.obsolete = True
            self.goNode.setObsolete(self.obsolete)

        elif name == "typedef":
            self.inTypedef = False

        elif name == "alt_id":
            self.synonyms[self.cdata.strip()] = self.goid

        elif name == "term" and not self.obsolete and self.namespace == self.elementNamespace:
            
            '''# not sure what the follow trying to do
            if not self.elementNamespace in self.namespaces.keys():
                self.namespaces[self.elementNamespace] = []

            self.namespaces[self.elementNamespace].append(self.goid)
            '''
            
            self.graph.add_node(self.goid, data=self.goNode)
            for parent in self.parents:
                '''
                currently the only edges being added are to the children of the current node
                '''
                if not parent in self.graph:
                    self.graph.add_node(parent, data=GONode(parent))
                self.graph.add_edge(parent, self.goid, relationship="parent_of")

        elif name == "name" and not self.obsolete and not self.goid == None and self.goNode.name == None: # and not self.name.has_key(self.goid):
            self.goNode.setName(self.cdata.strip())

        elif name == "defstr" and  not self.obsolete and not self.goid == None: # and not self.name.has_key(self.goid):
            self.goNode.setDescription(self.cdata.strip())

        elif name == "term" and self.obsolete and self.namespace == self.elementNamespace:
            self.obsoletes.add( self.goid )

    def characters(self, data):
        self.cdata += data

class GONode():
    '''from gographer.GONode.py'''
    def __init__ (self, goid=None, namespace=None, parents=None, obsolete=False,
                  name=None, description=None, genes=set(), propGenes=set(),
                  pmids=set(), propPmids=set(), wordVector=dict(), descendantCount=None,
                  mergedGenes=set(), mergedPmids=set(), mergedCount=0, infoLoss=0):
        self.goid = goid
        self.namespace = namespace
        self.parents = parents
        self.obsolete = obsolete
        self.name = name
        self.description = description
        self.genes = genes
        self.propGenes = propGenes
        self.pmids = pmids
        self.propPmids = propPmids
        self.wordVector = wordVector
        self.descendantCount = descendantCount
        self.mergedGenes = mergedGenes
        self.mergedPmids = mergedPmids
        self.mergedCount = mergedCount
        self.infoLoss = infoLoss


    ##Set the GO ID of the node
    # @param    goid    The GO ID that will be assigned to the node
    def setGOID (self, goid):
        self.goid = goid

    ##Returns the GO ID of the node
    def getGOID (self):
        return self.goid

    ##Set the namespace of the node
    # @param    namespace   The namespace that will be assigned to the node 
    def setNamespace (self, namespace):
        self.namespace = namespace

    ##Returns the namespace of the node
    def getNamespace (self):
        return self.namespace

    ##Set the parents of the node
    # @param    parents  The parents that will be assigned to the node
    def setParents (self, parents):
        self.parents = parents
        
    ##Returns the parents of the node
    def getParents (self):
        return self.parents

    ##Set the obsolete status of the node
    # @param    obsolete    The obsolete status that will be assigned to the node
    def setObsolete (self, obsolete):
        self.obsolete = obsolete

    ##Returns the obsolete status of the node (whether or not the term is obsolete)
    def getObsolete (self):
        return self.obsolete

    ##Set the name of the node
    # @param    name    The name that will be assigned to the node
    def setName (self, name):
        self.name = name

    ##Returns the name of the node
    def getName (self):
        return self.name
    
    ##Set the description of the node
    # @param    description The description that will be assigned to the node
    def setDescription (self, description):
        self.description = description

    ##Returns the description of the node
    def getDescription (self):
        return self.description

    ##Set the PubMed IDs associated with the node
    # @param    pmids   The set containing the PubMed IDs that are associated with the node.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPMIDs (self, pmids):
        self.pmids = pmids

    ##Returns the PubMed IDs associated with the node
    def getPMIDs (self):
        return self.pmids

    ##Adds list containing one or more PubMed ID tuples to the existing set of PubMed IDs
    # @param    pmid    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPMIDs(self, pmid):
        self.pmids = self.pmids.union(pmid)

    ##Set the propagated PubMed IDs associated with the node or its descendants
    # @param    propPmids   The set containing the propagated PubMed IDs that are associated with the node or its descendants.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPropagatedPMIDs (self, propPmids):
        self.propPmids = propPmids

    ##Returns the propagated PubMed IDs associated with the node or its descendants
    def getPropagatedPMIDs (self):
        return self.propPmids

    ##Adds list containing PubMed ID tuples to the existing set of propagated PubMed IDs
    # @param    pmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPropagatedPMIDs(self, pmids):
        self.propPmids = self.propPmids.union(pmids)

    ##Set the genes associated with the node
    # @param    genes   The set containing the genes that are associated with the node.
    #                   The genes are in tuples where it's the gene ID followed by the qualifier.
    def setGenes (self, genes):
        self.genes = genes

    ##Returns the genes associated with the node
    def getGenes (self):
        return self.genes

    ##Adds list containing one or more gene tuples to the existing set of genes
    # @param    genes    A list where each entry is a gene tuple, where it's the gene ID followed by the qualifier
    def addGenes(self, genes):
        self.genes = self.genes.union(genes)

    ##Set the propagated genes associated with the node or its descendants
    # @param    propGenes   The set containing the propagated PubMed IDs that are associated with the node or its descendants.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPropagatedGenes (self, propGenes):
        self.propGenes = propGenes

    ##Returns the propagated PubMed IDs associated with the node or its descendants
    def getPropagatedGenes (self):
        return self.propGenes

    ##Adds list containing PubMed ID tuples to the existing set of propagated PubMed IDs
    # @param    pmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPropagatedGenes(self, genes):
        self.propGenes = self.propGenes.union(genes)

    ##Calculates the word vector and stores this information
    # @param    corpus  The corpus that contains the information on the PubMed article
    # @param    pmids   A list of PubMed ID tuples to be added to the word vector, where it's the PubMed ID number followed by the qualifier
    #                   The propagated PubMed IDs will be used if none is given
    # @param    tokenizer   The tokenizer function that will be used on the text, a simple tokenizer is used if none is given.
    #                       Should take a string as an input, and outputs a string with words that are lower case and separated by a space
    # @param    stemmer The stemmer function that will be used to stem the words, the porter stemmer is used if none is given
    #                   Takes a word as an input and reports a stemmed word as an output.
    # @param    stopwords   A list of stop words that will not be included in the word vector, either as a list or a StopwordList.
    #                       An empty list is used if no stop word list is given.
    def calculateWordVector(self, corpus, tokenizer=Tokenizer().tokenize_word,
                            stemmer=PorterStemmer().stem, stopwords=[], pmids=None):
        wordVector = {}

        if not corpus:
            print "Not corpus is given, so a word vector can not be calculated."

        else:
            if not pmids:
                pmids = self.getPropagatedPMIDs()

            for pmid, qualifier in pmids:
                #Checks to make sure the PubMed article is in the corpus
                #and the qualifier does not contain the word 'NOT'
                if pmid in corpus.docs and "NOT" not in qualifier:
                    #Adds each word to the word vector
                    pmidWordVector = corpus.docs[pmid].getWordVector(tokenizer, stemmer, stopwords)
                    for word in pmidWordVector:
                        #if word not in stopwords:
                        if word in wordVector:
                            wordVector[word] += pmidWordVector[word]
                        else:
                            wordVector[word] = pmidWordVector[word]
        self.wordVector = wordVector
    
    ##Returns the word vector for the node calculated using propagated PubMed IDs
    # @param    corpus  The corpus that contains the information on the PubMed article
    def getWordVector(self, corpus=None, stopwords=[]):
        #Calculates the word vector if the current word vector is empty
        if len(self.wordVector) == 0:
            if len(self.getPropagatedPMIDs()) > 0:
                self.calculateWordVector(corpus, stopwords=stopwords)
        return self.wordVector

    ##Returns the number of descendants of the current node
    def getDescendantCount(self):
        return self.descendantCount

    
    ##Set the genes associated with the nodes that have been merged into the current node
    # @param    mergedGenes   The set containing the genes that are associated with the merged nodes.
    #                   The genes are in tuples where it's the gene ID followed by the qualifier.
    def setMergedGenes (self, mergedGenes):
        self.mergedGenes = mergedGenes

    ##Returns the genes associated with the nodes that have been merged into the current node
    def getMergedGenes (self):
        return self.mergedGenes

    ##Adds list containing one or more gene tuples to the existing set of merged genes
    # @param    mergedGenes    A list where each entry is a gene tuple, where it's the gene ID followed by the qualifier
    def addMergedGenes(self, mergedGenes):
        self.mergedGenes = self.mergedGenes.union(mergedGenes)


    ##Set the PubMed IDs associated with the nodes that have been merged into the current node
    # @param    mergedPmids   The set containing the PubMed IDs that are associated with the merged nodes.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setMergedPMIDs (self, mergedPmids):
        self.mergedPmids = mergedPmids

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getMergedPMIDs (self):
        return self.mergedPmids

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    mergedPmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addMergedPMIDs(self, mergedPmids):
        self.mergedPmids = self.mergedPmids.union(mergedPmids)

    ##Set the count of the number of nodes that have been merged into the current node
    # @param    mergedCount   The integer count of the number of nodes that have been merged into the current node
    def setMergedCount (self, mergedCount):
        self.mergedCount = mergedCount

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getMergedCount (self):
        return self.mergedCount

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    mergedPmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addMergedCount(self, mergedCount=1):
        self.mergedCount += mergedCount

    ##Set the PubMed IDs associated with the nodes that have been merged into the current node
    # @param    infoLoss   The set containing the PubMed IDs that are associated with the merged nodes.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setInfoLoss (self, infoLoss):
        self.infoLoss = infoLoss

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getInfoLoss (self):
        return self.infoLoss

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    infoLoss    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addInfoLoss(self, infoLoss):
        self.infoLoss += infoLoss
        
class GOGraph(networkx.DiGraph):
    '''from gographer.GOGraph.py'''
    ## Constructor.
    # @param    namespace     The branch of the GO ontology stored by this graph
    # @param    XMLFileName   The file contains the GO definition in
    #                         the format of the OBO in XML
    def __init__(self, namespace=None, GOOboXmlFileName=None):        
        DiGraph.__init__(self)

        self.namespace = namespace

        if GOOboXmlFileName != None:
            self.parseOboXml(GOOboXmlFileName)
        
    ## Parses the given OBO XML file and creates a GOGraph
    # @param    GOOboXMLFileName    The name of the OBO XML file to be parsed
    def parseOboXml(self, GOOboXmlFileName):
        parser = make_parser()
        handler = GOOboXmlHandler(self)
        parser.setContentHandler(handler)
        f = IOTools.openFile(GOOboXmlFileName, 'r')
        parser.parse(f)
        f.close()
        self.synonyms = handler.synonyms
        self.obsoletes = handler.obsoletes

    ## Returns the minimum depth of the node
    # @param    goid    The GO ID of the node whose depth will be found and returned
    def getLevel(self, goid):
        parents = self.predecessors(goid)
        if len(parents) == 0:
            return 0
        
        found = False
        level = 1
        while not found:
            nextDepth = list()
            for i in parents:
                pred = self.predecessors(i)
                if len(pred) == 0:
                    found = True
                    return level
                else:
                    nextDepth = list(set(nextDepth) | set(pred))
            parents = nextDepth
            level = level + 1


    ## Returns the namespace of the graph
    def getNameSpace(self):
        return self.namespace;

    ## Save self using pickle protocol.
    # @param filename The filename of the pickle to use
    def savePickle(self, filename="gograph.pickle"):
        try:
            f = open(filename, 'wb')
            cPickle.dump(self, f, protocol = cPickle.HIGHEST_PROTOCOL)
            f.close()
        except:
            print "Could not pickle graph"

    ## Load a pickle from the filesystem
    # @param    filename    The location of the pickle file to load
    @classmethod
    def loadPickle (klass, filename="gograph.pickle"):
        try:
            f = open(filename, 'rb')
            g = cPickle.load(f)
            return g
        except:
            print "Could not load pickle"
            return

    ## Get the description of a node
    # @param    goid   The GO ID of the node to get the description of
    def getNodeDescription(self, goid):
        if goid in self.nodes():
            return self.node[goid]['data'].getDescription()
        else:
            print 'Invalid goid'
            raise

    ## Set the description of a node
    # @param    goid   The GO ID of the node to get the description of
    # @param    descrip The description that will be assigned to the node
    def __setNodeDescription(self, goid, descrip):
        if goid in self.nodes():
            self.node[goid]['data'].setDescription(descrip)
        else:
            print 'Invalid goid'
            raise

    ## Get the description of a node
    # @param    goid   The GO ID of the node to get the description of
    def getNodeNamespace(self, goid):
        if goid in self.nodes():
            return self.node[goid]['data'].getNamespace()
        else:
            print 'Invalid goid'
            raise

    ## Set the description of a node
    # @param    goid   The GO ID of the node to get the description of
    # @param    namespace The namespace that will be assigned to the node
    def __setNodeNamespace(self, goid, namespace):
        if goid in self.nodes():
            self.node[goid]['namespace'] = namespace
        else:
            print 'Invalid goid'
            raise

    ## Calculates and stores the number of descendants each node in the graph has
    def calcDescendantCount(self):
        nodes = dict()
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()
        for node in sortedNodes:
            if node not in nodes:
                ids = set()
            else:
                ids = nodes[node]

            self.node[node]['data'].descendantCount = len(ids)
            ids = ids.union([node])

            for parent in self.predecessors(node):
                if parent not in nodes:
                    nodes[parent] = ids
                else:
                    nodes[parent] = nodes[parent].union(ids)

    ## Return the number of descendants of the given node. Calculates the count if it had not already been done so
    # @param    goid    The GO ID of the node to return the descendant count of
    def getDescendantCount(self, goid):
        if goid in self:
            count = self.node[goid]['data'].getDescendantCount()
            if not count:
                self.calcDescendantCount()
                return self.node[goid]['data'].getDescendantCount()
            else:
                return count



def getAncestors( g, node ):
    p = g.predecessors( node )
    ancestors = set(p)
    work = p[:]
    while work:
        x = work.pop()
        p = g.predecessors( x )
        work.extend( [ x for x in p if x not in ancestors ])
        ancestors.update(p)

    return ancestors

def readPValues( infile, synonyms = {} ):
    '''read pvalues and optionally fold changes from infile.'''
    term2pvalue, term2log2fold = {}, collections.defaultdict( float )

    def _parse2(data):
        goid = synonyms.get( data[0], data[0] )
        term2pvalue[goid] = float(data[1])

    def _parse3(data):
        goid = synonyms.get( data[0], data[0] )
        term2pvalue[goid] = float(data[1])
        term2log2fold[goid] = math.log( float(data[2]), 2 )

    def _parseGO(data):
        '''parse output from GO.py'''
        if len(data) < 14: return
        goid = synonyms.get( data[11], data[11] )
        term2pvalue[goid] = float( data[8] )
        term2log2fold[goid] = math.log( float(data[7])+0.0000001, 2 )

    parser = None
    for line in infile:
        if line.startswith("code"): 
            parser = _parseGO
            continue
        elif line.startswith("goid"):
            data = line[:-1].split("\t")
            if len(data) == 2:
                parser = _parse2
            else:
                parser = _parse3
            continue
        if parser == None:
            raise ValueError("don't know how to parse this file at line %s" % line)
        
        parser( line[:-1].split() )

    return term2pvalue, term2log2fold

def readSimrelMatrix( infile ):
    '''read simrel matrix from file.'''
    matrix, row_headers, col_headers = IOTools.readMatrix( infile )
    assert row_headers == col_headers
    return matrix, row_headers

def readTermFrequencies( infile ):
    '''read term frequencies and counts from file.'''
    counts, p = {}, {}
    for line in infile:
        if line.startswith("goid"):continue
        goid, f, c = line[:-1].split("\t")
        p[goid] = float(f)
        counts[goid] = set(c.split(";"))

    return counts, p

def readAncestors( infile ):
    '''read term frequencies and counts from file.'''
    ancestors = {}
    for line in infile:
        if line.startswith("goid"):continue
        goid, anc = line[:-1].split("\t")
        ancestors[goid] = set(anc.split(";"))
    return ancestors

def computeAncestors( graph, nodes ):
    '''build dictionary with node ancestors in graph.'''

    iteration, total = 0, len(nodes)

    E.info("computing ancestral information for %i nodes" % len(nodes))

    ancestors = {}
    for node in nodes:
        a = getAncestors( graph, node )
        # only take nodes with annotations
        a = set([ x for x in a if x in nodes ])
        # add self
        # a.add( node )
        ancestors[node] = a 
        iteration += 1
        if iteration % 500 == 0: 
            E.debug( "total = %i (%5.2f%%)" % (iteration, 100.0 * iteration / total) )
    return ancestors

def getRoot( graph ):
    '''return root of graph.'''
    root = graph.nodes()[0]
    p = graph.predecessors( root )
    while p != []:
        root = p[0]
        p = graph.predecessors( root )
    return root

def buildFilename( ontology, section, suffix ):
    return E.getOutputFile( "%s.%s.%s" % (ontology,section,suffix) )

def computeTermFrequencies( graph, go2genes ):

    nodes = graph.nodes()
    
    root = getRoot( graph )

    ############################################################
    # compute annotation counts - number of times a term
    # occurs in the background database
    E.info("computing annotation counts" )
    counts = collections.defaultdict( set )

    ############################################################
    # compute term frequencies p[term]
    # eq. 1 in Schlicker et al. (2006)
    E.info("computing term frequencies" )
    # compute counts
    for node in traversal.dfs_postorder( graph, source = root ):
    # older networkx:
    # for node in traversal.dfs_postorder_nodes( graph, source = root ):
        counts[node] = set(go2genes[node])
        for x in graph.successors( node ):
            counts[node].update( counts[x] )

    # normalize
    rf = float(len(counts[root]))
    p = {}
    for term, v in counts.iteritems():
        p[term] = len(v) / rf
        
    E.info( "computed %i frequencies" % len(p) )

    # debugging code - check that frequencies are increasing
    # monotonically.
    for node in nodes:
        pre = graph.predecessors( node )
        while pre != []:
            assert p[pre[0]] >= p[node]
            node = pre[0]
            pre = graph.predecessors( node )

    return counts, p

def buildSimrelMatrix( go2genes, go2info, ancestors, p ):
    '''build simrel matrix from an ontology DAG in *graph*
    and a mapping of go terms to genes.

    This method follows the procedure and syntax
    according to Schlicker et al. (2006) 
    BMC Bioinformatics 7:302
    '''

    ############################################################
    # Build symmetric similarity matrix between terms

    # only work on nodes with annotations
    lp = {}
    for term, v in p.iteritems():
        # ignore p = 1.0 as these will be downweighted to 
        # 0 (1.0 - p[c])
        if v > 0 and v < 1.0: lp[term] = math.log( v )

    test_nodes = lp.keys()
    node2index = dict( [ (y,x) for x,y in enumerate( test_nodes ) ] )

    # compute dictionary of ancestors for
    # each node. A simple dfs does not work as the
    # graph is not a tree.
    E.info("building %i x %i simrel matrix" % (len(test_nodes), len(test_nodes)))
    matrix = numpy.zeros( (len(test_nodes), len(test_nodes)) )
    iteration = 0

    total = len(test_nodes) * (len(test_nodes) - 1) / 2

    for c1, c2 in itertools.combinations( test_nodes, 2 ):
        iteration += 1
        if iteration % 1000000 == 0: 
            E.info( "total = %i (%5.2f%%)" % (iteration, 100.0 * iteration / total) )

        s1 = ancestors[c1]
        s2 = ancestors[c2]
        
        # print "c1=", c1, p[c1], lp[c1]
        # outputList( s1, go2info, p, lp )

        # print "c2=", c2, p[c2], lp[c2]
        # outputList( s2, go2info, p, lp )

        common = s1.intersection( s2 )

        i,j = node2index[c1], node2index[c2]

        f = 2.0 / (lp[c1] + lp[c2])
        # add pseudocount to avoid division by zero 
        vals = [ lp[c] * (1.0 - p[c] ) for c in common if c in lp ]
        if len(vals) > 0:
            # use min as vals are multiplied by negative factor
            val = f * min(vals)
            assert 0 <= val <= 1
            matrix[i][j] = val
            matrix[j][i] = val

    return matrix, test_nodes

def clusterSimrelMatrix( matrix, terms, ancestors, p, counts, term2pvalue, go2info, options):
    '''apply revigo clustering scheme to *matrix*.

    *ancestors*, *p*, and *counts* contain auxiliary data

    ancestors: dictionary mapping a term to all its ancestors
    p: dictionary mapping a term to its annotation frequency
    counts: dictionary mapping a term to all the genes it is annotated with
    go2info: dictionary mapping a term to its description
    term2pvalue: dictionary mapping a term to its p-value

    returns (matrix, terms)
    '''

    # cluster or load clustered simrel matrix
    E.info("working on simrel matrix with %i terms" % len(terms))
    term2index = dict( [ (y,x) for x,y in enumerate( terms ) ] )

    # shrink matrix to only those terms that are in result set
    terms = [ x for x in terms if x in term2pvalue ]

    E.info("after removing terms with no p-value: %i terms" % len(terms))
    if len(terms) == 0:
        E.warn( "no terms given after p-value filtering" )
        return None, [], []

    matrix = matrix[ numpy.array( [ term2index[x] for x in terms ]) ]
    matrix = matrix[ :, numpy.array( [ term2index[x] for x in terms ]) ]
    term2index = dict( [ (y,x) for x,y in enumerate( terms ) ] )

    # cluster according to revigo
    clusters = dict( [ (x, [x]) for x in terms ] )
    nrows = len(terms)

    counter = E.Counter()

    # row/column indices of members
    member_indices = []

    while 1:
        # find most similar pair of GO terms ti and tj with similarity simrel
        index = numpy.argmax( matrix )
        i = index // nrows
        j = index % nrows

        simrel = matrix[i][j]
        assert simrel == matrix.max()            

        # stop procedure if similarity is below threshold
        if simrel == 0:
            E.debug("similarity reached 0 - done" )
            break

        if simrel < options.max_simrel_threshold: 
            E.debug("similarity %f < %f - done" % (simrel, options.max_simrel_threshold))
            break

        nclusters = len(clusters)
        if nclusters <= options.min_clusters_threshold:
            E.debug("%i clusters reached at similarity %f - done" % (nclusters, simrel) )
            break

        t1 = terms[i]
        t2 = terms[j]
        E.debug( "clustering at similarity %f (%i,%i)" % (simrel, i,j))

        p1 = term2pvalue[t1]
        p2 = term2pvalue[t2]

        E.debug( "1. %s (%s): pvalue=%f, freq=%f" % (t1, go2info.get(t1, "na"), p1, p[t1] * 100.0))
        E.debug( "2. %s (%s): pvalue=%f, freq=%f" % (t2, go2info.get(t2, "na"), p2, p[t2] * 100.0))

        counter.tested += 1

        if p[t1] > options.min_frequency:
            # remove term with frequency > 5%
            counter.min_frequency += 1
            rep, mem = t2, t1

        elif p[t2] > options.min_frequency:
            # remove term with frequency > 5%
            counter.min_frequency += 1
            rep, mem = t1, t2

        elif term2pvalue[t1] < term2pvalue[t2]:
            # remove term with higher pvalue
            counter.min_pvalue += 1
            rep, mem = t1, t2

        elif term2pvalue[t2] < term2pvalue[t1]:
            # remove term with higher pvalue
            counter.min_pvalue += 1
            rep, mem = t2, t1

        elif t1 in ancestors[t2]:
            # t1 is parent of t2 - reject child term t2
            ovl = float(len(counts[t1].intersection( counts[t2] )))
            if (ovl / len(counts[t1])) >= options.min_child_threshold:
                # exception: if parent is > 75% child, keep child
                counter.ancestor_child += 1
                rep, mem = t2, t1
            else:
                counter.ancestor_parent += 1
                rep, mem = t1, t2

        elif t2 in ancestors[t1]:
            # t2 is parent of t1 - reject child term t1
            ovl = float(len(counts[t1].intersection( counts[t2] )))
            if (ovl / len(counts[t2])) >= options.min_child_threshold:
                # exception: if parent is > 75% child, keep child
                counter.ancestor_child += 1
                rep, mem = t1, t2
            else:
                counter.ancestor_parent += 1
                rep, mem = t2, t1

        else:
            # keep random term
            counter.random += 1
            rep, mem = min(t1,t2), max(t1,t2)

        # merge clusters
        clusters[rep].extend( clusters[mem] )
        del clusters[mem]
        # set values of mem to 0
        if mem == t1: 
            matrix[i,:] = 0
            matrix[:,i] = 0
            member_indices.append( i )
        else:
            matrix[j,:] = 0
            matrix[:,j] = 0
            member_indices.append( j )

    E.info( "%s" % str(counter))
    E.info( "after clustering: %i clusters" % len(clusters) )

    # shrink matrix
    member_indices = set( member_indices )
    rep_indices = [ x for x in xrange(nrows) if x not in member_indices ]
    terms = [ terms[x] for x in rep_indices ]
    rows = len(terms)
    matrix = matrix[ numpy.array( [ x for x in rep_indices ]) ]
    matrix = matrix[ :, numpy.array( [ x for x in rep_indices ]) ]

    return matrix, terms, clusters

@E.cachedfunction
def readGOGraph( filename_obo, namespace ):
    E.info( "reading GO graph from %s: namespace =%s" % (filename_obo, namespace))
    graph = GOGraph( namespace = namespace, GOOboXmlFileName=filename_obo )
    return graph

def main( argv ):

    parser = E.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"])

    parser.add_option("-o", "--filename-ontology", dest="filename_obo", type="string",
                      help="filename with ontology information in obo-xml format. "
                      " The latest version can always be"
                      " retrieved at `wget http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz`."
                      " [default=%default]." )

    parser.add_option("-g", "--filename-go", dest="filename_go", type="string",
                      help="filename with gene to go assignments "
                      " This is a tab-separated file with the columns <namespace> <gene_id> <goid> <description> <evidence>." 
                      " The evidence column is currently ignored. "
                      " [default=%default]." )

    parser.add_option("-p", "--filename-pvalues", dest="filename_pvalues", type="string",
                      help="filename with pvalues for each GO term "
                      " This is a tab-separated file with the columns <goid> <pvalue>."
                      " If a third column is present it will be interpreted as fold-change. "
                      " [default=%default]." )

    parser.add_option( "--ontology", dest="ontology", type="choice", action="append",
                       choices=("biol_process","cell_location","mol_function", "all" ),
                       help="ontologies to analyze. Ontologies are tested separately."
                       " 'all' will test all ontologies. "
                       " [default=%default]." )

    parser.add_option( "--max-similarity", dest="max_simrel_threshold", type="float",
                      help="cluster until simrel threshold is achieved [default=%default]." )

    parser.add_option( "--min-clusters", dest="min_clusters_threshold", type="float",
                      help="cluster only until at most # clusters remain [default=%default]." )

    parser.add_option( "--palette", dest="palette", type="choice",
                       choices=("rainbow", "gray", "blue-white-red",
                                "autumn", "bone", "cool", "copper", "flag", "gray", 
                                "hot", "hsv", "jet", "pink", "prism",
                                "spring", "summer", "winter", "spectral",
                                "RdBu", "RdGy", "BrBG", "BuGn", "Blues", "Greens", "Reds", "Oranges", "Greys" ),
                       help="colour palette [default=%default]")
    

    parser.add_option( "--reverse-palette", dest="reverse_palette", action="store_true",
                      help="reverse colour palette [default=%default]." )

    parser.set_defaults( filename_obo = 'go_daily-termdb.obo-xml.gz',
                         filename_go = "go.tsv.gz",
                         #filename_pvalues = "/ifs/projects/proj008/report/genelists.go.dir/rela_upstream-changed.biol_process.results", 
                         filename_pvalues = None,
                         filename_bg = None, 
                         filename_matrix = None,
                         filename_frequencies = None,
                         ontology = [],
                         # 0.9, 0.7, 0.5 and 0.4
                         # 0.4 very small clusters
                         max_simrel_threshold = 0.7,
                         min_frequency = 0.05,
                         min_pvalue = 0.05,
                         min_child_threshold = 0.75,
                         min_clusters_threshold = 0,
                         palette = "RdBu",
                         reverse_palette = False,
                         )

    options, args = E.Start( parser, argv, add_output_options = True )

    map_ontology2namespace = {
        'biol_process' :"biological_process",
        'cell_location': 'cellular_component',
        'mol_function': 'molecular_function' }

    if "all" in options.ontology:
        options.ontology = map_ontology2namespace.keys()

    #########################################
    if options.filename_pvalues == None:
        E.info("reading pvalues from stdin")
        all_term2pvalue, all_term2log2fold = readPValues( options.stdin )
    else:
        E.info("reading pvalues from %s " % options.filename_pvalues)
        with IOTools.openFile( options.filename_pvalues ) as infile:
            all_term2pvalue, all_term2log2fold = readPValues( infile )
        
    E.info( "read %i pvalue and %i fold assignments" % ( len(all_term2pvalue), len(all_term2log2fold) ))

    for test_ontology in options.ontology:
        
        E.info( "working on namespace %s" % test_ontology)

        namespace = map_ontology2namespace[test_ontology]

        ################################################################
        # compute or load ancestry information for all nodes
        fn = buildFilename( test_ontology, "ancestors", "tsv.gz" )
        if os.path.exists( fn ):
            E.info( "reading ancestor information from %s" % fn )
            with IOTools.openFile( fn ) as infile:
                ancestors = readAncestors( infile )
        else:
            graph = readGOGraph( options.filename_obo, namespace )
            ancestors = computeAncestors( graph, graph.nodes() )
            with IOTools.openFile( fn, "w" ) as outfile:
                outfile.write("goid\tancestors\n" )
                outfile.write("".join( [ "%s\t%s\n" % (x,";".join(v)) for x,v in ancestors.iteritems()] ) )

        #########################################
        graph = readGOGraph( options.filename_obo, namespace )

        #########################################
        E.info("reading gene to go assignments from %s " % options.filename_go)
        all_gene2gos, all_go2infos = GO.ReadGene2GOFromFile( IOTools.openFile( options.filename_go ), 
                                                             synonyms  = graph.synonyms,
                                                             obsolete = graph.obsoletes )

        #########################################
        # filter pvalues
        synonyms = graph.synonyms 
        term2pvalue = dict([ (synonyms.get(x,x),y) for x,y in all_term2pvalue.iteritems() ])
        term2log2fold = dict([ (synonyms.get(x,x),y) for x,y in all_term2log2fold.iteritems() ])

        if len(term2pvalue) == 0:
            E.warn( "no data - no output produced" )
            E.Stop()
            return
        
        go2info = all_go2infos[test_ontology]
        go2info = collections.defaultdict( str, dict( [(x.mGOId, x.mDescription) for x in go2info.values()] ) )
        gene2gos = all_gene2gos[test_ontology]

        if options.filename_bg:
            background = IOTools.readList( IOTools.openFile( filename_bg )) 
        else:
            background = list(gene2gos.keys())

        # invert mapping, at the same time insert missing ancestral counts
        go2genes = GO.buildGO2Genes( gene2gos, ancestors = ancestors )

        ################################################################
        # compute or load term frequences
        fn = buildFilename( test_ontology, "termfrequencies", "tsv.gz" )
        if os.path.exists( fn ):
            E.info( "reading term frequencies from %s" % fn )
            with IOTools.openFile( fn ) as infile:
                counts, p = readTermFrequencies( infile )
        else:
            counts, p = computeTermFrequencies( graph, go2genes )
            with IOTools.openFile( fn, "w" ) as outfile:
                outfile.write( "goid\tfrequency\tcounts\n" )
                for k in counts.keys():
                    outfile.write("%s\t%f\t%s\n" % (k, p[k], ";".join(counts[k]) ))

        # compute or load simrel matrix
        fn = buildFilename( test_ontology, "simrel", "matrix.gz" )
        if os.path.exists( fn ):
            E.info( "reading simrel matrix from %s" % fn )
            with IOTools.openFile( fn, "r") as infile:
                matrix, terms = readSimrelMatrix( infile )
        else:
            matrix, terms = buildSimrelMatrix( go2genes, go2info, ancestors, p )
            E.info("writing simrel matrix to file" )
            with IOTools.openFile( fn, "w" ) as outfile:
                IOTools.writeMatrix( outfile, matrix, terms, terms )
            
        fn = buildFilename( test_ontology, "cluster", "matrix.gz" )

        if os.path.exists( fn ):
            E.info( "reading clustered matrix matrix from %s" % fn )
            with IOTools.openFile( fn, "r") as infile:
                matrix, terms = readSimrelMatrix( infile )
        else:
            matrix, terms, clusters = clusterSimrelMatrix( matrix, terms, 
                                                           ancestors, p, counts, 
                                                           term2pvalue, go2info,
                                                           options )

            if matrix == None:
                E.warn( "%s: no output" % test_ontology)
                continue
            
            E.info( "after clustering: %i x %i matrix" % (len(terms), len(terms)))
            with IOTools.openFile( fn, "w" ) as outfile:
                IOTools.writeMatrix( outfile, matrix, terms, terms )

            fn = buildFilename( test_ontology, "cluster", "tsv")

            with IOTools.openFile( fn, "w" ) as outfile:
                outfile.write( "\t".join( ("rep", "mem", "pvalue", "frequency", "description") )+ "\n" )
                for rep, members in clusters.iteritems():
                    for mem in members:
                        outfile.write("%s\t%s\t%f\t%f\t%s\n" % (
                                rep, mem, term2pvalue[mem], p[mem], go2info.get(mem, "na")))


        # visualize output
        rows = len(terms)
            
        # revigo
        #     color = pvalue, size = frequency, egdewith = simrel score
        # me:
        #     size  = pvalue
        #     color = fold (blueish for depletion, reddish for enrichment), 
        #     edgewith = can't be done with networkx? (only edgecolor)

        # build graph from matrix
        term_graph = networkx.Graph()
        figure = plt.figure()

        term_graph.add_nodes_from( terms )
        for i,j in itertools.combinations( xrange( rows), 2 ):
            w = weight=matrix[i][j]
            if w < 0.1: continue
            term_graph.add_edge( terms[i], terms[j], weight = w )

        layout = networkx.spring_layout( term_graph )

        if options.reverse_palette:
            color_scheme = eval( "pylab.cm.%s_r" % options.palette)                    
        else:
            color_scheme = eval( "pylab.cm.%s" % options.palette)
        
        labels = dict( [ (x, go2info.get(x,x)) for x in terms] )
        node_size = [ -20.0 * math.log(term2pvalue[x],10) for x in term_graph.nodes() ]
        for x in term2pvalue.keys():
            term2log2fold[x] = numpy.random.normal(0, 5)
        node_color = [ term2log2fold[x] for x in term_graph.nodes() ]

        networkx.draw_networkx_edges( term_graph, layout, edge_color='0.50' )

        vmin = min(node_color)
        vmax = max(node_color)
        vmin = min( vmin, -vmax )
        vmax = max( vmax, -vmin )

        networkx.draw_networkx_nodes( term_graph, layout, 
                                      node_size = node_size,
                                      node_color = node_color,
                                      cmap = color_scheme,
                                      alpha = 0.8,
                                      vmin = vmin,
                                      vmax = vmax,
                                      linewidths = 0.5 )

        positions = dict( [ (x, numpy.array( (y[0],y[1]-0.02)) ) for x,y in layout.iteritems() ] )
        
        networkx.draw_networkx_labels( term_graph, 
                                       positions, 
                                       labels = labels,
                                       font_size = 5,
                                       alpha = 0.5 )

        # turn of tick marks
        gca = figure.gca()
        gca.set_xticklabels( [], visibile = False )
        gca.set_yticklabels( [], visibile = False )
        
        fn = buildFilename( test_ontology, "graph", "svg" )
        plt.colorbar()
        plt.savefig( fn )
        
        fn = buildFilename( test_ontology, "graph", "tsv.gz" )
        with IOTools.openFile( fn, "w") as outfile:
            outfile.write("goid\tx\ty\tpvalue\tl2fold\tdescription\n" )
            for term in terms:
                outfile.write("%s\t%f\t%f\t%f\t%f\t%s\n" % \
                                  (term, 
                                   layout[term][0],
                                   layout[term][1],
                                   term2pvalue[term],
                                   term2log2fold[term],
                                   go2info.get(term,"na") ) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
