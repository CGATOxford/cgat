################################################################################
#   Gene prediction pipeline 
#
#   $Id: GO.py 2883 2010-04-07 08:46:22Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, getopt, time, optparse, math, tempfile, subprocess, random
import collections

import scipy
import scipy.stats
import scipy.special
import numpy

import Database
import Experiment
import IOTools

USAGE="""program $Id: GO.py 2883 2010-04-07 08:46:22Z andreas $

calculate over-/under-representation of GO catergories in gene lists.
"""

# The following code was taken from:
#
# http://mail.python.org/pipermail/python-list/2006-January/359797.html
#
#
def lnchoose(n, m):
    nf = scipy.special.gammaln(n + 1)
    mf = scipy.special.gammaln(m + 1)
    nmmnf = scipy.special.gammaln(n - m + 1)
    return nf - (mf + nmmnf)

def hypergeometric_gamma(k, n1, n2, t):
    if t > n1 + n2:
        t = n1 + n2
    if k > n1 or k > t:
        return 0
    elif t > n2 and ((k + n2) < t):
        return 0
    else:
        c1 = lnchoose(n1,k)
        c2 = lnchoose(n2, t - k)
        c3 = lnchoose(n1 + n2 ,t)

    # print "hyperg:", k, n1, n2, t, math.exp(c1 + c2 - c3)
    return math.exp(c1 + c2 - c3)

def hypergeometric_P( k, n0, n1, t):

    GSL_DBL_EPSILON=1e-10

    assert t <= (n0+n1), "t larger than population size"
    assert n0 >= 0, "n0 < 0"
    assert n1 >= 0, "n1 < 0"

    if k >= n0 or k >= t:
        P = 1.0
    elif (k < 0.0):
        P = 0.0
    else:
        P = 0.0;
        mode = int( float(t*n0) / float(n0+n1))
        relerr = 1.0;
        if k < mode :
            i = k;
            relerr = 1.0;
            while(i >= 0 and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1
        else:
            i = mode
            relerr = 1.0;
            while(i <= k and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
            i = mode - 1;
            relerr = 1.0;
            while( i >=0 and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1
    return P

def hypergeometric_Q( k, n0, n1, t):

    GSL_DBL_EPSILON=1e-10

    assert t <= (n0+n1), "t larger than population size"    
    assert n0 >= 0, "n0 < 0"
    assert n1 >= 0, "n1 < 0"

    if k >= n0 or k >= t:
        P = 1.0
    elif (k < 0.0):
        P = 0.0
    else:
        P = 0.0;
        mode = int(float(t * n0) / float(n0+n1))
        relerr = 1.0
        if k < mode :
            i = mode
            relerr = 1.0
            while(i <= t and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
            i = mode - 1
            relerr = 1.0;
            while( i > k and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1
            
        else:
            i = k + 1
            relerr = 1.0
            while(i <= t and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma( i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
    return P


class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class GOEntry:
    
    mNameSpaceMap= {
        'molecular_function' : 'mol_function',
        'cellular_component' : 'cell_location',
        'biological_process' : 'biol_process', 
        }
    
    def __init__(self):
        pass
    
    def fromOBO( self, infile ):
        """read entry form an OBO formatted file."""

        self.mIsA = []

        while 1:
            line = infile.readline()
            if not line or line[0] == "\n": break
            data = line[:-1].split(":")
            term = data[0]
            rest = ":".join( data[1:] ).strip()
            if term == "name": self.mName = rest
            elif term == "id": self.mId = rest
            elif term == "namespace": self.mNameSpace = self.mNameSpaceMap[rest]
            elif term == "def": self.mDefinition = rest
            elif term == "exact_synonym": self.mSynonym = rest
            elif term == "is_a": self.mIsA.append( rest )
            elif term == "comment": self.mComment = rest
            elif term == "is_obsolete": self.mIsObsolete = True
            
def readOntology( infile ):
    """read ontology in OBO format from infile.

    returns a dictionary of Ontology entries.
    """
    result = {}
    while 1:
        line = infile.readline()
        if not line: break
        if line.startswith( "[Term]" ):
            go = GOEntry()
            go.fromOBO( infile )
            result[go.mId] = go
    
    return result

##-------------------------------------------------------------------------------
class GOSample:
    """store results from sampling.
    """
    def __init__(self, mmin, mmax, mmean, mstddev, mprobovers, mprobunders, counts):

        self.mMin = mmin
        self.mMax = mmax
        self.mMean = mmean
        self.mStddev = mstddev
        self.mProbabilitiesOverRepresentation = mprobovers
        self.mProbabilitiesUnderRepresentation = mprobunders
        self.mCounts = counts
##-------------------------------------------------------------------------------
class GOResult:

    mIsOverRepresented = False
    mGOId = None
    mSampleCountsCategory = 0
    mBackgroundCountsCategory = 0
    mSampleCountsTotal = 0
    mBackgroundCountsTotal = 0
    mProbabilityOverRepresentation = 0
    mProbabilityUnderRepresentation = 0

    def __init__(self, goid):
        self.mGOId = goid

    def UpdateProbabilities( self ):
        """calculate probabilities for given counts.

        """
        if self.mBackgroundCountsTotal == 0:
            return

        # various sanity checs
        assert self.mBackgroundCountsCategory >= self.mSampleCountsCategory, \
            "%s: more counts in foreground (%i) than in the background (%i) - make sure the foreground is part of the background." %\
            (self.mGOId, self.mSampleCountsCategory, self.mBackgroundCountsCategory )

        assert self.mBackgroundCountsTotal >= self.mBackgroundCountsCategory, \
            "%s: background: more counts in catagory (%i) than in total (%i)." %\
            (self.mGOId, self.mBackgroundCountsCategory, self.mBackgroundCountsTotal)

        assert self.mSampleCountsTotal >= self.mSampleCountsCategory, \
            "%s: forerground: more counts in catagory (%i) than in total (%i)." %\
            (self.mGOId, self.mSampleCountsCategory, self.mSampleCountsTotal)

        if self.mSampleCountsCategory == 0:
            self.mProbabilityOverRepresentation = 1.0
        else:
            self.mProbabilityOverRepresentation = hypergeometric_Q( self.mSampleCountsCategory - 1,
                                                                    self.mBackgroundCountsCategory,
                                                                    self.mBackgroundCountsTotal - self.mBackgroundCountsCategory,
                                                                    self.mSampleCountsTotal )
        
        self.mProbabilityUnderRepresentation = hypergeometric_P( self.mSampleCountsCategory,
                                                                 self.mBackgroundCountsCategory,
                                                                 self.mBackgroundCountsTotal - self.mBackgroundCountsCategory,
                                                                 self.mSampleCountsTotal )
        

        if self.mSampleCountsTotal == 0 or self.mBackgroundCountsCategory == 0:
            self.mRatio = "na"
        else:
            self.mRatio = float(self.mSampleCountsCategory) * self.mBackgroundCountsTotal / self.mSampleCountsTotal / self.mBackgroundCountsCategory 

    def __str__(self):
        """return string representation."""        
        return "%i\t%i\t%s\t%i\t%i\t%s\t%s\t%5.2e\t%5.2e" % (self.mSampleCountsCategory,
                                                             self.mSampleCountsTotal,
                                                             IOTools.prettyPercent( self.mSampleCountsCategory, self.mSampleCountsTotal ),
                                                             self.mBackgroundCountsCategory,
                                                             self.mBackgroundCountsTotal,
                                                             IOTools.prettyPercent( self.mBackgroundCountsCategory, self.mBackgroundCountsTotal ),
                                                             IOTools.prettyFloat( self.mRatio ),
                                                             self.mProbabilityOverRepresentation,
                                                             self.mProbabilityUnderRepresentation )

##-------------------------------------------------------------------------------                                           
class GOResults:

    def __init__(self):
        self.mResults = {}
        self.mNumGenes = 0
        self.mBackgroundCountsTotal = 0
        self.mSampleCountsTotal = 0

    def __str__(self):
        """return string representation."""
        lines = []
        lines.append( "\t".join(map(str, (self.mNumGenes, self.mBackgroundCountsTotal, self.mSampleCountsTotal))))
        for k, v in self.mResults.items():
            lines.append("%s\t%s" % (k, str(v)))
        return "\n".join(lines)

##-------------------------------------------------------------------------------
class GOInfo:
    mGOId = None
    mGOType = None
    mDescription = None
    
    def __init__( self,
                  goid = None,
                  go_type = None,
                  description = None ):

        self.mDescription = description
        self.mGOId = goid
        self.mGOType = go_type

    def __str__(self):
        return "\t".join(map(str, (self.mGOId, self.mGOType, self.mDescription)))
    

##-------------------------------------------------------------------------------
class GOMatch(GOInfo):
    mEvidence = None
    
    def __init__( self,
                  goid = None,
                  go_type = None,
                  description = None,
                  evidence = None ):

        GOInfo.__init__(self, goid, go_type, description)
        self.mEvidence = evidence

    def __str__(self):
        return "\t".join(map(str, (self.mGOId, self.mGOType, self.mDescription, self.mEvidence)))

##---------------------------------------------------------------------
def FilterByGOIds(  gene2go, go2info ):

    """
    filter gene_id to go_id lookup by a list of go_ids

    returns a new gene2go mapping.
    
    used to restrict GO terms to GO_slim and remove alternates

    gene2go 				# starting set, map of genes to go terms
    go2info				# alt ids are repeats of superceding ids
    """

    filtered_gene2go = {}

    for gene_id in gene2go.keys():
        new_go = set()
        for go in gene2go[gene_id]:
            if go.mGOId in go2info:
                new_go.add( go )

        if new_go:
            filtered_gene2go[gene_id] = list(new_go)

    return filtered_gene2go

##---------------------------------------------------------------------
def MapGO2Slims(  gene2go, go2slim, ontology = None ):
    """filter gene2go lookup by a list of go_ids in go2slim.

    gene2go: map of genes to go terms
    go2slim: map of go categories to goslim go categories   

    If ontology is given, missing descriptions of go entries
    are added from the ontology.

    returns a new gene2go mapping.
    """

    ## build map of go identifiers to go info
    map_go2info = {}
    if ontology:
        for go in ontology.values():
            map_go2info[go.mId] = GOInfo( goid=go.mId, 
                                          go_type = go.mNameSpace,
                                          description = go.mName )
    else:
        for gene_id, gos in gene2go.items():
            for go in gos:
                map_go2info[go.mGOId] = go

    filtered_gene2go = {}

    for gene_id, gos in gene2go.items():    
        new_go = set()
        for go in gos:
            if go.mGOId in go2slim:
                for gg in go2slim[go.mGOId]:
                    if gg in map_go2info:
                        new_go.add( map_go2info[gg] )
                    else:
                        raise IndexError, "description for mapped go term not present: %s -> %s" % (go.mGOId, gg)
        if new_go:
            filtered_gene2go[gene_id] = list(new_go)

    return filtered_gene2go

##------------------------------------------------------------------------------
def GetGOSlims( infile ):
    """
    returns a map of go identifiers to slim categories

    Input is the output of Chris Mungal's map2slim.pl.
    """

    go2go = {}
    for line in infile:
        if line[:len("part_of")] == "part_of": continue

        mapped, parents = line.split("//")
        go, goslims = mapped.split("=>")
        goslims = goslims.split(" ")
        if len(goslims) == 0 : continue

        go2go[go.strip()] = filter( lambda x: len(x), map(lambda x: x.strip(), goslims))
        
    return go2go

##------------------------------------------------------------------------------
def GetGOFrequencies( gene2go, genes ):
    """count number of each go category in gene list.

    return a tuple containing:
    * the total number of GO categories found.
    * dictionary of counts per GO category
    * dictionary of genes found with GO categories
    """
    counts = {}
    total = 0
    found_genes = {}
    
    for gene_id in genes:
        
        if gene_id not in gene2go: continue

        found_genes[gene_id] = 1
        for go in gene2go[gene_id]:
            if go.mGOId not in counts: counts[go.mGOId] = 0
            counts[go.mGOId] += 1
            total += 1

    return total, counts, found_genes

##------------------------------------------------------------------------------
def AnalyseGO( gene2go,
               genes,
               genes_background = (),
               do_probabilities = True):
    """analyse go ids.

    goids: list of goids to analyse
    genes: sample set of genes
    genes_background: background set of genes (default: all)
    """
    
    if genes_background == ():
        genes_background = gene2go.keys()

    result = GOResults()
    
    ## get background frequencies
    (background_counts_total, background_counts, background_genes) = \
        GetGOFrequencies( gene2go, 
                          genes_background )
    
    result.mBackgroundCountsTotal = background_counts_total
    result.mBackgroundNumCategories = len(background_counts)
    result.mBackgroundGenes = background_genes

    ## get sample frequencies
    (sample_counts_total, sample_counts, sample_genes) = \
                          GetGOFrequencies( gene2go, 
                                            genes )

    result.mNumGenes = len(genes)

    result.mSampleCountsTotal = sample_counts_total
    result.mSampleNumCategories = len(sample_counts)
    result.mSampleGenes = sample_genes
    
    # test for over or underrepresented categories in the slims
    # report results for all go categories in the background
    # so that also categories completely absent in the foreground (sample)
    # are considered.
    for go_id in background_counts.keys():
        
        result_go = GOResult(go_id)

        if go_id in sample_counts:
            # use gene counts
            result_go.mSampleCountsCategory = sample_counts[go_id]
        else:
            result_go.mSampleCountsCategory = 0

        result_go.mSampleCountsTotal = len(sample_genes)
        result_go.mBackgroundCountsTotal = len(background_genes)
        result_go.mBackgroundCountsCategory = background_counts[go_id]

        if do_probabilities:
            try:
                result_go.UpdateProbabilities()
            except AssertionError, msg:
                print msg
                print "# error while calculating probabilities for %s" % go_id
                print "# genes in sample", sample_genes
                print "# counts in sample: %i / %i" % ( result_go.mSampleCountsCategory, result_go.mSampleCountsTotal)
                print "# counts in background %i / %i" % (result_go.mBackgroundCountsCategory, result_go.mBackgroundCountsTotal)
                for x in sample_genes.keys():
                    for y in gene2go[x]:
                        print x, str(y)

                sys.exit(0)
        
        result.mResults[go_id] = result_go
        
    return result

##--------------------------------------------------------------------------------
def GetGOStatement( go_type, database, species ):
    """build statement to get GO assignments for genes from ENSEMBL."""
    
    if database in ("ensembl_mart_27_1", ):
        statement = """SELECT DISTINCTROW
        gene_stable_id, glook_%s_id, description, olook_evidence_code
        FROM %s.%s_gene_ensembl__go_%s__look
        WHERE glook_%s_id IS NOT NULL
        GROUP BY gene_stable_id, glook_%s_id, description
        ORDER BY gene_stable_id
        """ % (go_type,
               database, species, go_type,
               go_type, go_type)
        
    elif database in ("ensembl_mart_31", "ensembl_mart_37", "ensembl_mart_41" ):
        statement = """SELECT DISTINCTROW
        gene_stable_id, glook_%s_id, description, olook_evidence_code
        FROM %s.%s_go_%s__go_%s__main
        WHERE glook_%s_id IS NOT NULL
        GROUP BY gene_stable_id, glook_%s_id, description
        ORDER BY gene_stable_id
        """ % (go_type,
               database, species, go_type, go_type,
               go_type, go_type)

    elif re.search( "core", database):

        if go_type == "biol_process":
            go_type = "biological_process"
        elif go_type == "mol_function":
            go_type = "molecular_function"
        elif go_type == "cell_location":
            go_type = "cellular_component"
        else:
            raise "unknown go_type %s" % go_type

        x = re.search("(\d+)", database)
        if not x:
            raise "can't find version number in database %s" % database

        version = int(x.groups()[0])
        if version <= 54:
            go_database = "ensembl_go_%s" % version
            go_field = "acc"
            statement = """SELECT DISTINCTROW
        g.stable_id, xref.dbprimary_acc, go.name, 'NA'
        FROM gene, transcript, translation, 
        gene_stable_id as g, object_xref as o, xref,
        %(go_database)s.term AS go
        WHERE gene.gene_id = transcript.gene_id
        AND transcript.transcript_id = translation.transcript_id
        AND g.gene_id = gene.gene_id
        AND translation.translation_id = o.ensembl_id
        AND xref.xref_id = o.xref_id
        AND go.%(go_field)s = xref.dbprimary_acc
        AND go.term_type = '%(go_type)s'
        AND xref.external_db_id = 1000
        """ % locals()

        else:
            go_database = "ensembl_ontology_%s" % version
            go_field = "accession"

            statement = """SELECT DISTINCTROW
        g.stable_id, xref.dbprimary_acc, go.name, 'NA'
        FROM gene, transcript, translation, 
        gene_stable_id as g, object_xref as o, xref,
        %(go_database)s.term AS go,
        %(go_database)s.ontology AS ontology
        WHERE gene.gene_id = transcript.gene_id
        AND transcript.transcript_id = translation.transcript_id
        AND g.gene_id = gene.gene_id
        AND translation.translation_id = o.ensembl_id
        AND xref.xref_id = o.xref_id
        AND go.%(go_field)s = xref.dbprimary_acc
        AND go.ontology_id = ontology.ontology_id 
        AND ontology.namespace = '%(go_type)s'
        AND xref.external_db_id = 1000
        """ % locals()

    else:
        raise "unknown ensmart version %s" % database

    return statement

##--------------------------------------------------------------------------------
def ReadGene2GOFromDatabase( dbhandle, go_type, database, species ):
    """read go assignments from ensembl database.

    returns a dictionary of lists.
    (one to many mapping of genes to GO categories)
    and a dictionary of go-term to go information

    Note: assumes that external_db_id for GO is 1000
    """

    statement = GetGOStatement( go_type, database, species )
    result = dbhandle.Execute(statement).fetchall()

    gene2go = {}
    go2info = {}
    for gene_id, goid, description, evidence in result:
        gm = GOMatch( goid, go_type, description, evidence )
        gi = GOInfo( goid, go_type, description )
        if gene_id not in gene2go: gene2go[gene_id] = []
        gene2go[gene_id].append(gm)
        go2info[goid] = gi
    
    return gene2go, go2info

##---------------------------------------------------------------------------
def DumpGOFromDatabase( outfile, 
                        dbhandle, 
                        options ):

    """read go assignments from database.

    and dump them into a flatfile.
    (one to many mapping of genes to GO categories)
    and a dictionary of go-term to go information
    """

    if options.loglevel >= 1:    
        options.stdlog.write("# category\ttotal\tgenes\tcategories\n" )

    all_genes = collections.defaultdict( int )
    all_categories = collections.defaultdict( int )
    all_ntotal = 0

    outfile.write("go_type\tgene_id\tgo_id\tdescription\tevidence\n" )

    for go_type in options.go_category:

        genes = collections.defaultdict( int )
        categories = collections.defaultdict( int )
        ntotal = 0

        statement = GetGOStatement( go_type, options.database, options.species )
        
        results = dbhandle.Execute(statement).fetchall()

        for result in results:
            outfile.write( "\t".join(map( str, (go_type,) + result))+ "\n")
            gene_id, goid, description, evidence = result
            genes[gene_id] += 1
            categories[goid] += 1
            ntotal += 1
            all_genes[gene_id] += 1
            all_categories[goid] += 1
            all_ntotal += 1
            
        if options.loglevel >= 1:    
            options.stdlog.write( "# %s\t%i\t%i\t%i\n" % (go_type, ntotal,
                                                          len(genes),
                                                          len(categories) ) )


    if options.loglevel >= 1:    
        options.stdlog.write( "# %s\t%i\t%i\t%i\n" % ("all", 
                                                      all_ntotal,
                                                      len(all_genes),
                                                      len(all_categories) ) )

    return 

##---------------------------------------------------------------------------
def ReadGene2GOFromFile( infile ):
    """reads GO mappings for all go_types from a
    file.

    returns two maps: gene2go maps genes to go categories
    and go2info maps go categories to information.
    """

    gene2gos = {}
    go2infos = {}
    for line in infile:
        if line[0] == "#": continue
        go_type, gene_id, goid, description, evidence = line[:-1].split("\t")
        if go_type == "go_type": continue

        gm = GOMatch( goid, go_type, description, evidence )
        gi = GOInfo( goid, go_type, description )
        if go_type not in gene2gos:
            gene2gos[go_type] = {}
            go2infos[go_type] = {}

        gene2go = gene2gos[go_type]
        go2info = go2infos[go_type]
        
        if gene_id not in gene2go: gene2go[gene_id] = []
        gene2go[gene_id].append(gm)
        go2info[goid] = gi
    
    return gene2gos, go2infos

##---------------------------------------------------------------------------
def CountGO( gene2go ):
    """count number of genes and go categories in mapping."""

    cats = {}
    nmaps = 0
    for k,vv in gene2go.items():
        for v in vv:
            nmaps += 1                    
            cats[v.mGOId] = 1
            
    return len(gene2go), len(cats), nmaps

##---------------------------------------------------------------------------
def countGOs( gene2gos ):
    """return map of number of genes and go categories in mapping."""
    genes, goids = collections.defaultdict( int ), collections.defaultdict( int )

    for cat, gene2go in gene2gos.iteritems():
        for gene_id, vv in gene2go.iteritems():
            genes[gene_id] += 1
            for v in vv:
                goids[v.mGOId] += 1
    return genes, goids

##---------------------------------------------------------------------------
def ReadGeneList( filename_genes, options ):
    
    """read gene list from filename."""

    if filename_genes == "-":
        infile = sys.stdin
    else:
        infile = open(filename_genes,"r")

    genes = map( lambda x: x[:-1].split("\t")[0], filter( lambda x: x[0] != "#", infile.readlines()))
    infile.close()
    if options.loglevel >= 1:        
        print "# read %i genes from %s" % (len(genes), filename_genes)

    ## apply transformation
    if options.gene_pattern:
        rx = re.compile(options.gene_pattern)
        genes = map( lambda x: rx.search( x ).groups()[0], genes )
            
    #############################################################
    ## make non-redundant
    xx = {}
    for x in genes: xx[x] = 1
    genes = xx.keys()

    if filename_genes != "-":
        infile.close()

    if options.loglevel >= 1:
        print "# after filtering: %i nonredundant genes." % (len(genes))

    return genes

##---------------------------------------------------------------------------
def GetCode( v ):
    """return a code for over/underrepresentation."""

    if v.mProbabilityOverRepresentation < v.mProbabilityUnderRepresentation:
        code = "+"
    elif v.mProbabilityOverRepresentation > v.mProbabilityUnderRepresentation:
        code = "-"
    else:
        code = "?"
    return code

##---------------------------------------------------------------------------    
def convertGo2Goslim( options ):
    """read gene list with GO assignments and convert to GO slim categories."""
    
    Experiment.info( "reading GO assignments from stdin" )
    gene2gos, go2infos = ReadGene2GOFromFile( options.stdin )
    input_genes, input_goids = countGOs( gene2gos )

    #############################################################
    ## read GO ontology from file
    assert options.filename_ontology, "please supply a GO ontology"
    Experiment.info( "reading ontology from %s" % (options.filename_ontology) )
        
    infile = open(options.filename_ontology)
    ontology = readOntology( infile )
    infile.close()
        
    go2infos = collections.defaultdict( dict )
    # substitute go2infos
    for go in ontology.values():
        go2infos[go.mNameSpace][go.mId] = GOInfo( go.mId,
                                                  go_type = go.mNameSpace,
                                                  description = go.mName )

    Experiment.info( "reading GO assignments from %s" % options.filename_slims)
    go_slims = GetGOSlims( open(options.filename_slims, "r") )

    if options.loglevel >=1:
        v = set()
        for x in go_slims.values():
            for xx in x: v.add(xx)
        Experiment.info( "read go slims from %s: go=%i, slim=%i" %\
                                  ( options.filename_slims,
                                    len(go_slims), 
                                    len( v ) ))
                
    output_goids, output_genes = set(), set()
    noutput = 0
    options.stdout.write( "\t".join( ("go_type", "gene_id", "go_id", "description", "evidence" ) ) + "\n" )
    for category, gene2go in gene2gos.iteritems():
        gene2go = MapGO2Slims( gene2go, go_slims, ontology )
        for gene_id, values in gene2go.iteritems():
            output_genes.add( gene_id )
            for go in values:
                output_goids.add( go.mGOId)
                options.stdout.write( "%s\t%s\t%s\t%s\t%s\n" % \
                                          ( go.mGOType,
                                            gene_id,
                                            go.mGOId,
                                            go.mDescription,
                                            "NA", ) )
                noutput += 1

    Experiment.info( "ninput_genes=%i, ninput_goids=%i, noutput_gene=%i, noutput_goids=%i, noutput=%i" % \
                         (len(input_genes), len(input_goids),
                          len(output_genes), len(output_goids),
                          noutput) )
##---------------------------------------------------------------------------    
if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: GO.py 2883 2010-04-07 08:46:22Z andreas $", usage=USAGE)

    dbhandle = Database.Database()
    
    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use [default=%default]." )

    parser.add_option("-i", "--slims", dest="filename_slims", type="string",
                      help="filename with GO SLIM categories [default=%default].")

    parser.add_option( "-g", "--genes", dest="filename_genes", type="string",
                       help="filename with genes to analyse [default=%default]." )

    parser.add_option( "-b", "--background", dest="filename_background", type="string",
                       help="filename with background genes to analyse [default=%default]." )

    parser.add_option( "-o", "--sort-order", dest="sort_order", type="choice",
                       choices=("fdr", "pover", "ratio" ),
                       help="output sort order [default=%default]." )

    parser.add_option( "-c", "--category", dest="go_category", type="string",
                       help="go category to analyse [biol_process|cell_location|mol_function] [default=%default]." )

    parser.add_option( "-t", "--threshold", dest="threshold", type="float",
                       help="significance threshold [>1.0 = all ]. If --fdr is set, this refers to the fdr, otherwise it is a cutoff for p-values." )

    parser.add_option ("--filename-dump", dest="filename_dump", type="string",
                       help="dump GO category assignments into a flatfile [default=%default]." )

    parser.add_option ("--filename-ontology", dest="filename_ontology", type="string",
                       help="filename with ontology in OBO format [default=%default]." )

    parser.add_option ( "--filename-input", dest="filename_input", type="string",
                       help="read GO category assignments from a flatfile [default=%default]." )

    parser.add_option ( "--sample", dest="sample", type="int",
                       help="do sampling (with # samples) [default=%default]." )

    parser.add_option ( "--filename-output-pattern", dest = "output_filename_pattern", type="string",
                        help="pattern with output filename pattern (should contain: %(go)s and %(section)s ) [default=%default]")

    parser.add_option ( "--output-filename-pattern", dest = "output_filename_pattern", type="string",
                        help="pattern with output filename pattern (should contain: %(go)s and %(section)s ) [default=%default]")
    
    parser.add_option ( "--fdr", dest="fdr", action="store_true",
                       help="calculate and filter by FDR [default=%default]." )

    parser.add_option ( "--go2goslim", dest="go2goslim", action="store_true",
                       help="convert go assignments in STDIN to goslim assignments and write to STDOUT [default=%default]." )

    parser.add_option ( "--gene-pattern", dest = "gene_pattern", type="string",
                        help="pattern to transform identifiers to GO gene names [default=%default].")

    parser.add_option( "--filename-map-slims", dest="filename_map_slims", type="string",
                       help="write mapping between GO categories and GOSlims [default=%default].")

    parser.add_option( "--get-genes", dest="get_genes", type="string",
                       help="list all genes in the with a certain GOID [default=%default]." )
    
    parser.set_defaults( species = None,
                         filename_genes = "-",
                         filename_background = None,
                         filename_slims = None,
                         go_category = "biol_process,cell_location,mol_function",
                         filename_categories = None,
                         filename_dump = None,
                         sample = 0,
                         fdr = False,
                         output_filename_pattern = None,
                         threshold = 0.05,
                         filename_map_slims = None,
                         gene_pattern = None,
                         sort_order = "ratio",
                         get_genes = None )

    (options, args) = Experiment.Start( parser, add_mysql_options = True )
    options.go_category = options.go_category.split(",")

    if options.fdr and options.sample == 0:
        print USAGE
        raise "please supply a sample size for determining the empirical FDR"

    if options.go2goslim:
        convertGo2Goslim( options )
        Experiment.Stop()
        sys.exit(0)

    #############################################################
    ## dump GO
    if options.filename_dump:
        if options.loglevel >= 1:        
            options.stdlog.write( "# dumping GO categories to %s\n" % (options.filename_dump) )
            sys.stdout.flush()

        dbhandle.Connect( options )
            
        outfile = open(options.filename_dump, "w")
        DumpGOFromDatabase( outfile,
                            dbhandle,
                            options )
        outfile.close()
        Experiment.Stop()
        sys.exit(0)

    #############################################################
    ## read GO categories from file
    if options.filename_input:
        if options.loglevel >= 1:        
            options.stdlog.write( "# reading association of categories and genes from %s\n" % (options.filename_input) )
            sys.stdout.flush()
        
        infile = open(options.filename_input)
        gene2gos, go2infos = ReadGene2GOFromFile( infile )
        infile.close()

    #############################################################
    ## read GO ontology from file
    if options.filename_ontology:
        if options.loglevel >= 1:        
            options.stdlog.write( "# reading ontology from %s\n" % (options.filename_ontology) )
            sys.stdout.flush()
        
        infile = open(options.filename_ontology)
        ontology = readOntology( infile )
        infile.close()
        
        go2infos = collections.defaultdict( dict )
        ## substitute go2infos
        for go in ontology.values():
            go2infos[go.mNameSpace][go.mId] = GOInfo( go.mId,
                                                      go_type = go.mNameSpace,
                                                      description = go.mName )

    #############################################################
    ## get foreground gene list
    genes = ReadGeneList( options.filename_genes, options )
        
    #############################################################
    ## get background
    if options.filename_background:
        background = ReadGeneList( options.filename_background, options )
        missing = set(genes).difference( set(background))
        assert len(missing) == 0, "%i genes in foreground but not in background: %s" % (len(missing), str(missing))
    else:
        background = ()

    #############################################################
    ## get go categories for genes
    for go_category in options.go_category:

        #############################################################
        ## get/read association of GO categories to genes
        if options.filename_input:
            gene2go, go2info = gene2gos[go_category], go2infos[go_category]
        else:
            if options.loglevel >= 1:
                options.stdlog.write( "# reading data from database ..." )
                sys.stdout.flush()

            dbhandle.Connect( options )
            gene2go, go2info = ReadGene2GOFromDatabase( dbhandle,
                                                        go_category,
                                                        options.database, options.species )

            if options.loglevel >= 1:
                options.stdlog.write( "finished.\n" )
                sys.stdout.flush()

        if options.loglevel >= 1:
            ngenes, ncategories, nmaps = CountGO( gene2go )        
            options.stdlog.write( "# read GO assignments: %i genes mapped to %i categories (%i maps)\n" % (ngenes, ncategories, nmaps) )

        #############################################################
        ## sanity check:            
        ## are all of the foreground genes in the dataset
        ## missing = set(genes).difference( set(gene2go.keys()) )
        ## assert len(missing) == 0, "%i genes in foreground set without GO annotation: %s" % (len(missing), str(missing))

        #############################################################            
        ## read GO slims and map GO categories to GO slim categories
        if options.filename_slims:
            go_slims = GetGOSlims( open(options.filename_slims, "r") )
            
            if options.loglevel >=1:
                v = set()
                for x in go_slims.values():
                    for xx in x: v.add(xx)
                options.stdlog.write( "# read go slims from %s: go=%i, slim=%i\n" %\
                                          ( options.filename_slims,
                                            len(go_slims), 
                                            len( v ) ))

                                       

            if options.filename_map_slims:
                if options.filename_map_slims == "-":
                    outfile = options.stdout
                else:
                    outfile=open(options.filename_map_slims, "w")

                outfile.write( "GO\tGOSlim\n" )
                for go, go_slim in go_slims.items():
                    outfile.write("%s\t%s\n" % (go, go_slim))

                if outfile != options.stdout:
                    outfile.close()

            gene2go = MapGO2Slims( gene2go, go_slims, ontology )

            if options.loglevel >=1:
                ngenes, ncategories, nmaps = CountGO( gene2go )
                options.stdlog.write( "# after go slim filtering: %i genes mapped to %i categories (%i maps)\n" % (ngenes, ncategories, nmaps) )

        #############################################################
        ## Just dump out the gene list
        if options.get_genes:
            fg, bg, ng = [], [], []

            for gene, vv in gene2go.items():
                for v in vv:
                    if v.mGOId == options.get_genes:
                        if gene in genes:
                            fg.append( gene )
                        elif gene in background:
                            bg.append( gene )
                        else:
                            ng.append( gene )

            ## skip to next GO class
            if not (bg or ng): continue
                
            options.stdout.write( "# genes in GO category %s\n" % options.get_genes )
            options.stdout.write( "gene\tset\n" )
            for x in fg: options.stdout.write("%s\t%s\n" % ("fg", x))
            for x in bg: options.stdout.write("%s\t%s\n" % ("bg", x))           
            for x in ng: options.stdout.write("%s\t%s\n" % ("ng", x))                       

            if options.loglevel >= 1:
                options.stdlog.write("# nfg=%i, nbg=%i, nng=%i\n" % (len(fg), len(bg), len(ng) ))
                
            Experiment.Stop()
            sys.exit(0)
                  
        #############################################################################
        ## sampling
        ## for each GO-category:
        ##      get maximum and minimum counts in x samples -> calculate minimum/maximum significance
        ##      get average and stdev counts in x samples -> calculate z-scores for test set
        samples = {}
        simulation_min_pvalues = None
        if options.sample:

            sample_size = options.sample
            # List of all minimum probabilities in simulation
            simulation_min_pvalues = []
            if options.loglevel >= 1:
                options.stdlog.write( "# sampling: calculating %i samples: " % (sample_size))
                sys.stdout.flush()
                
            counts = {}
            prob_overs = {}
            prob_unders = {}
            
            for x in range(sample_size):

                if options.loglevel >= 1:
                    options.stdlog.write( "." )
                    options.stdlog.flush()
                    
                ## get shuffled array of genes from background
                sample_genes = random.sample( background, len(genes) )

                go_results = AnalyseGO( gene2go , sample_genes, background )

                pairs = go_results.mResults.items()
                
                for k, v in pairs:
                    if k not in counts:
                        counts[k] = []
                        prob_overs[k] = []
                        prob_unders[k] = []
                    
                    counts[k].append( v.mSampleCountsCategory )
                    prob_overs[k].append( v.mProbabilityOverRepresentation )
                    prob_unders[k].append( v.mProbabilityUnderRepresentation )                    
                    
                    simulation_min_pvalues.append( min( v.mProbabilityUnderRepresentation,
                                                        v.mProbabilityOverRepresentation ) )


            if options.loglevel >= 1:
                sys.stdout.write("\n")
                sys.stdout.flush()

            if options.loglevel >= 1:
                options.stdlog.write( "# sampling: sorting %i P-Values\n" % len(simulation_min_pvalues) )
                sys.stdout.flush()
            
            simulation_min_pvalues.sort()
            simulation_min_pvalues = numpy.array(simulation_min_pvalues)
                
            samples = {}

            if options.output_filename_pattern:

                filename = options.output_filename_pattern % { 'go': go_category, 'section': "samples" }
                options.stdlog.write( "# sampling results go to %s\n" % filename )
                outfile = open(filename, "w")
            else:
                outfile = sys.stdout
                
            outfile.write( "\t".join( ("goid", "min", "max", "mean", "median", "stddev", 
                                       "CI95lower", "CI95upper",
                                       "pover", "punder", "goid",
                                       "category", "description") ) + "\n" )
            for k in counts.keys():

                c = counts[k]

                prob_overs[k].sort()
                prob_unders[k].sort()

                s = GOSample(min(c),
                             max(c),
                             scipy.mean(c),
                             numpy.std(c),
                             numpy.array(prob_overs[k]),
                             numpy.array(prob_unders[k]),
                             counts[k] )
                
                samples[k] = s
                              
                if k in go2info:
                    n = go2info[k]
                else:
                    n = "?"
                
                outfile.write( "%s\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" %\
                               (k,
                                min(c),
                                max(c),
                                scipy.mean(c),
                                scipy.median(c),
                                numpy.std(c),
                                scipy.stats.scoreatpercentile( c, 5 ),
                                scipy.stats.scoreatpercentile( c, 95 ),
                                min(prob_overs[k]),
                                min(prob_unders[k]),
                                n ))
            if options.output_filename_pattern:
                outfile.close()
                

        #############################################################
        ## do the analysis
        go_results = AnalyseGO( gene2go, genes, background )

        pairs = go_results.mResults.items()

        #############################################################
        ## calculate fdr for each hypothesis
        fdrs = {}

        if options.fdr:

            if options.loglevel >= 1:
                options.stdlog.write( "# calculating the FDRs\n" )
                sys.stdout.flush()
                
            observed_min_pvalues = map( lambda x: min(x[1].mProbabilityOverRepresentation,
                                                      x[1].mProbabilityUnderRepresentation),
                                        pairs )
            observed_min_pvalues.sort()
            observed_min_pvalues = numpy.array( observed_min_pvalues )

            for k, v in pairs:

                if k in samples:
                    s = samples[k]
                else:
                    raise "category %s not in samples" % k
                
                ## calculate values for z-score
                if s.mStddev > 0:
                    zscore = abs(float(v.mSampleCountsCategory) - s.mMean) / s.mStddev
                else:
                    zscore = 0.0

                #############################################################
                # FDR:
                # For each p-Value p at node n:
                #   a = average number of nodes in each simulation run with P-Value < p
                #           this can be obtained from the array of all p-values and all nodes
                #           simply divided by the number of samples.
                #      aka: expfpos=experimental false positive rate
                #   b = number of nodes in observed data, that have a P-Value of less than p.
                #      aka: pos=positives in observed data
                #   fdr = a/b
                pvalue = min(v.mProbabilityOverRepresentation,
                             v.mProbabilityUnderRepresentation)

                # calculate values for FDR: 
                # nfdr = number of entries with P-Value better than node.
                a = 0
                while a < len(simulation_min_pvalues) and \
                          simulation_min_pvalues[a] < pvalue:
                    a += 1
                a = float(a) / float(sample_size)
                b = 0
                while b < len(observed_min_pvalues) and \
                        observed_min_pvalues[b] < pvalue:
                    b += 1

                if b > 0:
                    fdr = min(1.0, float(a) / float(b))
                else:
                    fdr = 0.0
                    
                fdrs[k] = (fdr, a, b)

        if options.sort_order == "fdr":
            pairs.sort( lambda x, y: cmp(fdrs[x[0]], fdrs[y[0]] ) )           
        elif options.sort_order == "ratio":
            pairs.sort( lambda x, y: cmp(x[1].mRatio, y[1].mRatio))
        elif options.sort_order == "pover":
            pairs.sort( lambda x, y: cmp(x[1].mProbabilityOverRepresentation, y[1].mProbabilityOverRepresentation))

        #############################################################
        ## output the result selected
        if options.output_filename_pattern:
            filename = options.output_filename_pattern % { 'go': go_category, 'section': "results" }
            if options.loglevel >= 1:
                options.stdlog.write( "# results go to %s\n" % filename)
            outfile = open(filename, "w")
        else:
            outfile = sys.stdout
        
        headers = ["code", "goid", "scount", "stotal", "spercent", "bcount", "btotal", "bpercent",
                   "ratio",
                   "pover", "punder", "goid", "category", "description"]

        if options.sample:
            headers += ["min", "max", "zscore", "mpover", "mpunder", "pos", "expfpos", "fdr", "nfdr_expected",
                        "expected", "CI95lower", "CI95upper" ]
        
        outfile.write("\t".join(headers) + "\n" )

        nselected = 0

        for k, v in pairs:
            
            code = GetCode( v )

            is_ok = True

            if options.fdr:
                (fdr, expfpos, pos) = fdrs[k]
                is_ok = fdr < options.threshold
            else:
                is_ok = min(v.mProbabilityOverRepresentation, v.mProbabilityUnderRepresentation) < options.threshold
                    
            if is_ok:
                
                nselected += 1
                
                if k in go2info:
                    n = go2info[k]
                else:
                    n = "?"
                    
                outfile.write("%s\t%s\t%s\t%s" % (code, k, str(v), n))

                if options.sample:

                    if k in samples:
                        s = samples[k]
                    else:
                        outfile.write("\n")
                    
                    ## calculate values for z-score
                    if s.mStddev > 0:
                        zscore = abs(float(v.mSampleCountsCategory) - s.mMean) / s.mStddev
                    else:
                        zscore = 0.0

                    # the number of expected false positives is the current FDR times the
                    # number of hypothesis selected.
                    nexpected = nselected * fdr
                    
                    outfile.write("\t%i\t%i\t%f\t%5.2e\t%5.2e\t%i\t%5.2f\t%5.2e\t%6.4f\t%6.4f\t%6.4f\t%6.4f" %\
                                  (s.mMin,
                                   s.mMax,
                                   zscore,
                                   min(s.mProbabilitiesOverRepresentation),
                                   min(s.mProbabilitiesUnderRepresentation),
                                   pos, expfpos, fdr, nexpected,
                                   scipy.mean( s.mCounts ),
                                   scipy.stats.scoreatpercentile( s.mCounts, 5 ),
                                   scipy.stats.scoreatpercentile( s.mCounts, 95 ),
                                   ) )
                outfile.write("\n")

        if options.output_filename_pattern:
            outfile.close()

        #############################################################
        ## output the full result
            
        if options.output_filename_pattern:
            filename = options.output_filename_pattern % { 'go': go_category, 'section': "overall" }
            if options.loglevel >= 1:
                options.stdlog.write( "# a list of all categories and pvalues goes to %s\n" % filename )
            outfile = open(filename, "w")
        else:
            outfile = sys.stdout
        
        headers = ["code", "goid", "passed", 
                   "scount", "stotal", "spercent", 
                   "bcount", "btotal", "bpercent",
                   "ratio", "pover", "punder", 
                   "goid", "category", "description"]

        if options.sample:
            headers += ["min", "max", "zscore", "mpover", "mpunder", "pos", "expfpos", "fdr", "mean", "CI95lower", "CI95upper" ]

        outfile.write("\t".join(headers) + "\n" )

        for k, v in pairs:
            
            code = GetCode( v )

            is_ok = "0"

            if options.fdr:
                (fdr, expfpos, pos) = fdrs[k]
                if fdr < options.threshold: is_ok = "1"
            else:
                if min(v.mProbabilityOverRepresentation, v.mProbabilityUnderRepresentation) < options.threshold: is_ok = "1"
                
            if k in go2info:
                n = go2info[k]
            else:
                n = "?"

            outfile.write("%s\t%s\t%s\t%s\t%s" % (code, k, is_ok, str(v), n))

            if options.sample:

                if k in samples:
                    s = samples[k]
                else:
                    outfile.write("\n")

                if options.fdr:
                    (fdr, expfpos, pos) = fdrs[k]
                else:
                    (fdr, expfpos, pos) = 0, 0, 0
                    
                ## calculate values for z-score
                if s.mStddev > 0:
                    zscore = abs(float(v.mSampleCountsCategory) - s.mMean) / s.mStddev
                else:
                    zscore = 0.0

                # the number of expected false positives is the current FDR times the
                # number of hypothesis selected.
                outfile.write("\t%i\t%i\t%f\t%5.2e\t%5.2e\t%i\t%5.2f\t%5.2e\t%f\t%f\t%f" %\
                                  (s.mMin,
                                   s.mMax,
                                   zscore,
                                   min(s.mProbabilitiesOverRepresentation),
                                   min(s.mProbabilitiesUnderRepresentation),
                                   pos, expfpos, fdr,
                                   scipy.mean( s.mCounts ),
                                   scipy.stats.scoreatpercentile( s.mCounts, 5 ),
                                   scipy.stats.scoreatpercentile( s.mCounts, 95 ),
                                   ) )
            outfile.write("\n")

        if options.output_filename_pattern:
            outfile.close()

        #############################################################
        ## output parameters
        ngenes, ncategories, nmaps = CountGO( gene2go )

        if options.output_filename_pattern:
            filename = options.output_filename_pattern % { 'go': go_category, 'section': "parameters" }
            if options.loglevel >= 1:
                options.stdlog.write( "# parameters go to %s\n" % filename )
            outfile = open(filename, "w")
        else:
            outfile = sys.stdout
            
        outfile.write( "# input go mappings for category '%s'\n" % go_category )
        outfile.write( "value\tparameter\n" )
        outfile.write( "%i\tmapped genes\n" % ngenes )
        outfile.write( "%i\tmapped categories\n" % ncategories )
        outfile.write( "%i\tmappings\n" % nmaps )

        nbackground = len(background)
        if nbackground == 0:
            nbackground = len(go_results.mBackgroundGenes)
            
        outfile.write( "%i\tgenes in sample\n" % len(genes) )
        outfile.write( "%i\tgenes in sample with GO assignments\n" % (len(go_results.mSampleGenes)) )
        outfile.write( "%i\tinput background\n" % nbackground )
        outfile.write( "%i\tgenes in background with GO assignments\n" % (len(go_results.mBackgroundGenes)) )
        outfile.write( "%i\tassociations in sample\n"     % go_results.mSampleCountsTotal )
        outfile.write( "%i\tassociations in background\n" % go_results.mBackgroundCountsTotal )
        outfile.write( "%s\tpercent genes in sample with GO assignments\n" % (IOTools.prettyPercent( len(go_results.mSampleGenes) , len(genes), "%5.2f" )))
        outfile.write( "%s\tpercent genes background with GO assignments\n" % (IOTools.prettyPercent( len(go_results.mBackgroundGenes), nbackground, "%5.2f" )))

        outfile.write( "%i\tsignificant results reported\n" % nselected )
        outfile.write( "%6.4f\tsignificance threshold\n" % options.threshold )        

        if options.output_filename_pattern:
            outfile.close()

        #############################################################
        ## output the fg patterns
            
        #############################################################
        ## Compute reverse map
        go2genes = {}

        for gene, gos in gene2go.items():
            if gene not in genes: continue
            for go in gos:
                if go.mGOId not in go2genes:
                    go2genes[go.mGOId] = []
                go2genes[go.mGOId].append( gene )
            
        if options.output_filename_pattern:
            filename = options.output_filename_pattern % { 'go': go_category, 'section': "fg" }
            if options.loglevel >= 1:
                options.stdlog.write( "# results go to %s\n" % filename )
            outfile = open(filename, "w")
        else:
            outfile = sys.stdout

        headers = ["code", "goid", "scount", "stotal", "spercent", "bcount", "btotal", "bpercent",
                   "ratio", 
                   "pover", "punder", "goid", "category", "description", "fg"]

        for k, v in pairs:

            code = GetCode( v )            

            if k in go2info:
                n = go2info[k]
            else:
                n = "?"
                
            if k in go2genes:
                g = ";".join( go2genes[k] )
            else:
                g = ""

            outfile.write("%s\t%s\t%s\t%s\t%s\n" % (code, k, str(v), n, g ) )

        if outfile != sys.stdout:
            outfile.close()

    Experiment.Stop()
