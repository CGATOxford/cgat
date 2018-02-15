'''
GO.py - compute GO enrichment from gene lists
=============================================

:Tags: Python

Code
----

'''
import sys
import re
import math
import random
import collections

import scipy
import scipy.stats
import scipy.special
import numpy
from CGAT import Stats as Stats
from CGAT import Experiment as E
from CGAT import IOTools as IOTools
from CGAT import Database as Database
from CGAT import CSV as CSV

from rpy2.robjects import r as R

MIN_FLOAT = sys.float_info.min


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
        c1 = lnchoose(n1, k)
        c2 = lnchoose(n2, t - k)
        c3 = lnchoose(n1 + n2, t)

    # print "hyperg:", k, n1, n2, t, math.exp(c1 + c2 - c3)
    return max(math.exp(c1 + c2 - c3), MIN_FLOAT)


def hypergeometric_P(k, n0, n1, t):

    GSL_DBL_EPSILON = 1e-10

    assert t <= (n0 + n1), "t larger than population size"
    assert n0 >= 0, "n0 < 0"
    assert n1 >= 0, "n1 < 0"

    if k >= n0 or k >= t:
        P = 1.0
    elif (k < 0.0):
        P = 0.0
    else:
        P = 0.0
        mode = int(float(t * n0) / float(n0 + n1))
        relerr = 1.0
        if k < mode:
            i = k
            relerr = 1.0
            while(i >= 0 and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1
        else:
            i = mode
            relerr = 1.0
            while(i <= k and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
            i = mode - 1
            relerr = 1.0
            while(i >= 0 and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1
    return P


def hypergeometric_Q(k, n0, n1, t):

    GSL_DBL_EPSILON = 1e-10
    assert t <= (n0 + n1), "t larger than population size"
    assert n0 >= 0, "n0 < 0"
    assert n1 >= 0, "n1 < 0"
    if k >= n0 or k >= t:
        P = 1.0
    elif (k < 0.0):
        P = 0.0
    else:
        P = 0.0
        mode = int(float(t * n0) / float(n0 + n1))
        relerr = 1.0
        if k < mode:
            i = mode
            relerr = 1.0
            while(i <= t and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
            i = mode - 1
            relerr = 1.0
            while(i > k and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i -= 1

        else:
            i = k + 1
            relerr = 1.0
            while(i <= t and relerr > GSL_DBL_EPSILON and P < 1.0):
                tmp = hypergeometric_gamma(i, n0, n1, t)
                P += tmp
                relerr = tmp / P
                i += 1
    return P


class Error(Exception):

    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)

    def _get_message(self, message):
        return self._message

    def _set_message(self, message):
        self._message = message
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

    mNameSpaceMap = {
        'molecular_function': 'mol_function',
        'cellular_component': 'cell_location',
        'biological_process': 'biol_process',
    }

    def __init__(self, default_namespace="ontology"):
        self.mNameSpace = default_namespace

    def fromOBO(self, section):
        """read entry form an OBO formatted file."""

        self.mIsA = []

        for line in section:

            data = line[:-1].split(":")
            term = data[0]
            rest = ":".join(data[1:]).strip()
            if term == "name":
                self.mName = rest
            elif term == "id":
                self.mId = rest
            elif term == "namespace":
                self.mNameSpace = self.mNameSpaceMap.get(rest, rest)
            elif term == "def":
                self.mDefinition = rest
            elif term == "exact_synonym":
                self.mSynonym = rest
            elif term == "is_a":
                self.mIsA.append(rest)
            elif term == "comment":
                self.mComment = rest
            elif term == "is_obsolete":
                self.mIsObsolete = True


def readOntology(infile):
    """read ontology in OBO format from infile.

    returns a dictionary of Ontology entries.
    """
    result = {}

    def iterate_blocks(infile):

        lines = []

        for line in infile:
            if line.strip() == "":
                if lines:
                    yield lines
                lines = []
                continue

            lines.append(line)

    default_namespace = "ontology"

    for section in iterate_blocks(infile):

        if section[0].startswith("[Term]"):
            go = GOEntry(default_namespace=default_namespace)
            go.fromOBO(section)
            result[go.mId] = go
        else:
            for line in section:
                data = line[:-1].split(":")
                if data[0] == "default-namespace":
                    default_namespace = data[1].strip()

    return result


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


class GOResult:

    mIsOverRepresented = False
    mGOId = None
    mSampleCountsCategory = 0
    mBackgroundCountsCategory = 0
    mSampleCountsTotal = 0
    mBackgroundCountsTotal = 0
    mProbabilityOverRepresentation = 0
    mProbabilityUnderRepresentation = 0
    mPValue = 1.0

    def __init__(self, goid=None):
        self.mGOId = goid

    def UpdateProbabilities(self):
        """calculate probabilities for given counts.

        """
        if self.mBackgroundCountsTotal == 0:
            return

        # various sanity checs
        assert self.mBackgroundCountsCategory >= self.mSampleCountsCategory, \
            "%s: more counts in foreground (%i) than in the background (%i) - make sure the foreground is part of the background." %\
            (self.mGOId, self.mSampleCountsCategory,
             self.mBackgroundCountsCategory)

        assert self.mBackgroundCountsTotal >= self.mBackgroundCountsCategory, \
            "%s: background: more counts in category (%i) than in total (%i)." %\
            (self.mGOId, self.mBackgroundCountsCategory,
             self.mBackgroundCountsTotal)

        assert self.mSampleCountsTotal >= self.mSampleCountsCategory, \
            "%s: forerground: more counts in category (%i) than in total (%i)." %\
            (self.mGOId, self.mSampleCountsCategory, self.mSampleCountsTotal)

        if self.mSampleCountsCategory == 0:
            self.mProbabilityOverRepresentation = 1.0
        else:
            self.mProbabilityOverRepresentation = hypergeometric_Q(self.mSampleCountsCategory - 1,
                                                                   self.mBackgroundCountsCategory,
                                                                   self.mBackgroundCountsTotal -
                                                                   self.mBackgroundCountsCategory,
                                                                   self.mSampleCountsTotal)

        self.mProbabilityUnderRepresentation = hypergeometric_P(self.mSampleCountsCategory,
                                                                self.mBackgroundCountsCategory,
                                                                self.mBackgroundCountsTotal -
                                                                self.mBackgroundCountsCategory,
                                                                self.mSampleCountsTotal)

        self.mPValue = min(
            self.mProbabilityOverRepresentation, self.mProbabilityUnderRepresentation)

        if self.mSampleCountsTotal == 0 or self.mBackgroundCountsCategory == 0:
            self.mRatio = "na"
        else:
            self.mRatio = float(self.mSampleCountsCategory) * self.mBackgroundCountsTotal / \
                self.mSampleCountsTotal / self.mBackgroundCountsCategory

    def getHeaders(self):
        return ["scount", "stotal", "spercent",
                "bcount", "btotal", "bpercent",
                "ratio",
                "pvalue", "pover", "punder"]

    def __str__(self):
        """return string representation."""
        return "%i\t%i\t%s\t%i\t%i\t%s\t%s\t%6.4e\t%6.4e\t%6.4e" % \
            (self.mSampleCountsCategory,
             self.mSampleCountsTotal,
             IOTools.prettyPercent(
                 self.mSampleCountsCategory, self.mSampleCountsTotal),
             self.mBackgroundCountsCategory,
             self.mBackgroundCountsTotal,
             IOTools.prettyPercent(
                 self.mBackgroundCountsCategory, self.mBackgroundCountsTotal),
             IOTools.val2str(self.mRatio),
             self.mPValue,
             self.mProbabilityOverRepresentation,
             self.mProbabilityUnderRepresentation)


class GOResults:

    '''container for go results.'''

    def __init__(self):
        # dictionary of (GOID,GoResult) tuples
        self.mResults = {}
        self.mNumGenes = 0
        self.mBackgroundCountsTotal = 0
        self.mSampleCountsTotal = 0

    def __str__(self):
        """return string representation."""
        lines = []
        lines.append("\t".join(
            map(str, (self.mNumGenes, self.mBackgroundCountsTotal, self.mSampleCountsTotal))))
        for k, v in list(self.mResults.items()):
            lines.append("%s\t%s" % (k, str(v)))
        return "\n".join(lines)


class GOInfo:
    mGOId = None
    mGOType = None
    mDescription = None

    def __init__(self,
                 goid=None,
                 go_type=None,
                 description=None):

        self.mDescription = description
        self.mGOId = goid
        self.mGOType = go_type

    def __str__(self):
        if self.mGOId is None:
            return "\t".join(map(str, ("", "", "")))
        else:
            return "\t".join(map(str, (self.mGOId, self.mGOType, self.mDescription)))

    def getHeaders(self):
        return ["goid", "go_catagory", "go_description"]


class GOMatch(GOInfo):
    mEvidence = None

    def __init__(self,
                 goid=None,
                 go_type=None,
                 description=None,
                 evidence=None):

        GOInfo.__init__(self, goid, go_type, description)
        self.mEvidence = evidence

    def __str__(self):
        return "\t".join(map(str, (self.mGOId, self.mGOType, self.mDescription, self.mEvidence)))


def FilterByGOIds(gene2go, go2info):
    """
    filter gene_id to go_id lookup by a list of go_ids

    returns a new gene2go mapping.

    used to restrict GO terms to GO_slim and remove alternates

    gene2go 				# starting set, map of genes to go terms
    go2info				# alt ids are repeats of superceding ids
    """

    filtered_gene2go = {}

    for gene_id in list(gene2go.keys()):
        new_go = set()
        for go in gene2go[gene_id]:
            if go.mGOId in go2info:
                new_go.add(go)

        if new_go:
            filtered_gene2go[gene_id] = list(new_go)

    return filtered_gene2go


def MapGO2Slims(gene2go, go2slim, ontology=None):
    """filter gene2go lookup by a list of go_ids in go2slim.

    gene2go: map of genes to go terms
    go2slim: map of go categories to goslim go categories

    If ontology is given, missing descriptions of go entries
    are added from the ontology.

    returns a new gene2go mapping.
    """

    # build map of go identifiers to go info
    map_go2info = {}
    if ontology:
        for go in list(ontology.values()):
            map_go2info[go.mId] = GOInfo(goid=go.mId,
                                         go_type=go.mNameSpace,
                                         description=go.mName)
    else:
        for gene_id, gos in list(gene2go.items()):
            for go in gos:
                map_go2info[go.mGOId] = go

    filtered_gene2go = {}

    for gene_id, gos in list(gene2go.items()):
        new_go = set()
        for go in gos:
            if go.mGOId in go2slim:
                for gg in go2slim[go.mGOId]:
                    if gg in map_go2info:
                        new_go.add(map_go2info[gg])
                    else:
                        raise IndexError(
                            "description for mapped go term not present: %s -> %s" %
                            (go.mGOId, gg))
        if new_go:
            filtered_gene2go[gene_id] = list(new_go)

    return filtered_gene2go


def GetGOSlims(infile):
    """
    returns a map of go identifiers to slim categories

    Input is the output of Chris Mungal's map2slim.pl.
    """

    go2go = {}
    for line in infile:
        if line[:len("part_of")] == "part_of":
            continue

        mapped, parents = line.split("//")
        go, goslims = mapped.split("=>")
        goslims = goslims.split(" ")
        if len(goslims) == 0:
            continue

        go2go[go.strip()] = [x for x in [x.strip() for x in goslims] if len(x)]

    return go2go


def GetGOFrequencies(gene2go, genes):
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

        if gene_id not in gene2go:
            continue

        found_genes[gene_id] = 1
        for go in gene2go[gene_id]:
            if go.mGOId not in counts:
                counts[go.mGOId] = 0
            counts[go.mGOId] += 1
            total += 1

    return total, counts, found_genes


def AnalyseGO(gene2go,
              genes,
              genes_background=None,
              do_probabilities=True):
    """analyse go ids.

    goids: list of goids to analyse
    genes: sample set of genes
    genes_background: background set of genes (default: all)
    """
    if genes_background is None:
        genes_background = list(gene2go.keys())

    result = GOResults()

    # get background frequencies
    (background_counts_total, background_counts, background_genes) = \
        GetGOFrequencies(gene2go,
                         genes_background)

    result.mBackgroundCountsTotal = background_counts_total
    result.mBackgroundNumCategories = len(background_counts)
    result.mBackgroundGenes = background_genes

    # get sample frequencies
    (sample_counts_total, sample_counts, sample_genes) = \
        GetGOFrequencies(gene2go,
                         genes)

    result.mNumGenes = len(genes)

    result.mSampleCountsTotal = sample_counts_total
    result.mSampleNumCategories = len(sample_counts)
    result.mSampleGenes = sample_genes

    # test for over or underrepresented categories in the slims
    # report results for all go categories in the background
    # so that also categories completely absent in the foreground (sample)
    # are considered.
    for go_id in list(background_counts.keys()):

        result_go = GOResult(go_id)

        # use gene counts
        result_go.mSampleCountsCategory = sample_counts.get(go_id, 0)
        result_go.mSampleCountsTotal = len(sample_genes)
        result_go.mBackgroundCountsTotal = len(background_genes)
        result_go.mBackgroundCountsCategory = background_counts[go_id]

        E.debug("processing %s: genes in foreground=%i, genes in backgound=%i, sample_counts=%i, background_counts=%i" %
                (go_id,
                 len(sample_genes),
                 len(background_genes),
                 sample_counts.get(go_id, 0),
                 background_counts.get(go_id, 0),
                 )
                )

        if do_probabilities:
            try:
                result_go.UpdateProbabilities()
            except AssertionError as msg:
                print(msg)
                print("# error while calculating probabilities for %s" % go_id)
                print("# genes in sample", sample_genes)
                print("# counts in sample: %i out of %i total" % (
                    result_go.mSampleCountsCategory, result_go.mSampleCountsTotal))
                print("# counts in background %i out of %i total" % (
                    result_go.mBackgroundCountsCategory, result_go.mBackgroundCountsTotal))
                for x in list(sample_genes.keys()):
                    for y in gene2go[x]:
                        print(x, str(y))

                sys.exit(0)

        result.mResults[go_id] = result_go

    return result


def GetGOStatement(go_type, database, species):
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

    elif database in ("ensembl_mart_31", "ensembl_mart_37", "ensembl_mart_41"):
        statement = """SELECT DISTINCTROW
        gene_stable_id, glook_%s_id, description, olook_evidence_code
        FROM %s.%s_go_%s__go_%s__main
        WHERE glook_%s_id IS NOT NULL
        GROUP BY gene_stable_id, glook_%s_id, description
        ORDER BY gene_stable_id
        """ % (go_type,
               database, species, go_type, go_type,
               go_type, go_type)

    elif re.search("core", database):

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

        elif version <= 66:
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

        elif version <= 88:
            go_database = "ensembl_ontology_%s" % version
            go_field = "accession"

            statement = """SELECT DISTINCTROW
        gene.stable_id, xref.dbprimary_acc, go.name, 'NA'
        FROM gene, transcript, translation,
        object_xref as o, xref,
        %(go_database)s.term AS go,
        %(go_database)s.ontology AS ontology
        WHERE gene.gene_id = transcript.gene_id
        AND transcript.transcript_id = translation.transcript_id
        AND translation.translation_id = o.ensembl_id
        AND xref.xref_id = o.xref_id
        AND go.%(go_field)s = xref.dbprimary_acc
        AND go.ontology_id = ontology.ontology_id
        AND ontology.namespace = '%(go_type)s'
        AND xref.external_db_id = 1000
        """ % locals()

        else:
            go_database = "ensembl_ontology_%s" % version
            go_field = "accession"

            statement = """SELECT DISTINCTROW
        gene.stable_id, xref.dbprimary_acc, go.name, 'NA'
        FROM gene, transcript,
        object_xref as o, xref,
        %(go_database)s.term AS go,
        %(go_database)s.ontology AS ontology
        WHERE gene.gene_id = transcript.gene_id
        AND transcript.transcript_id = o.ensembl_id
        AND o.ensembl_object_type = 'Transcript'
        AND xref.xref_id = o.xref_id
        AND go.%(go_field)s = xref.dbprimary_acc
        AND go.ontology_id = ontology.ontology_id
        AND ontology.namespace = '%(go_type)s'
        AND xref.external_db_id = 1000
        """ % locals()

    else:
        raise "unknown ensmart version %s" % database

    return statement


def ReadGene2GOFromDatabase(dbhandle, go_type, database, species):
    """read go assignments from ensembl database.

    returns a dictionary of lists.
    (one to many mapping of genes to GO categories)
    and a dictionary of go-term to go information

    Note: assumes that external_db_id for GO is 1000
    """

    statement = GetGOStatement(go_type, database, species)
    result = Database.executewait(dbhandle, statement,
                                  retries=0).fetchall()

    gene2go = {}
    go2info = collections.defaultdict(GOInfo)
    for gene_id, goid, description, evidence in result:
        gm = GOMatch(goid, go_type, description, evidence)
        gi = GOInfo(goid, go_type, description)
        if gene_id not in gene2go:
            gene2go[gene_id] = []
        gene2go[gene_id].append(gm)
        go2info[goid] = gi

    return gene2go, go2info


def DumpGOFromDatabase(outfile,
                       dbhandle,
                       options):
    """read go assignments from database.

    and dump them into a flatfile.
    (one to many mapping of genes to GO categories)
    and a dictionary of go-term to go information
    """

    E.info("category\ttotal\tgenes\tcategories")

    all_genes = collections.defaultdict(int)
    all_categories = collections.defaultdict(int)
    all_ntotal = 0

    outfile.write("go_type\tgene_id\tgo_id\tdescription\tevidence\n")

    for go_type in options.ontology:

        genes = collections.defaultdict(int)
        categories = collections.defaultdict(int)
        ntotal = 0
        statement = GetGOStatement(go_type, options.database_name,
                                   options.species)

        results = Database.executewait(
            dbhandle, statement, retries=0).fetchall()

        for result in results:
            outfile.write("\t".join(map(str, (go_type,) + result)) + "\n")
            gene_id, goid, description, evidence = result
            genes[gene_id] += 1
            categories[goid] += 1
            ntotal += 1
            all_genes[gene_id] += 1
            all_categories[goid] += 1
            all_ntotal += 1

        E.info("%s\t%i\t%i\t%i" % (go_type, ntotal,
                                   len(genes),
                                   len(categories)))

    E.info("%s\t%i\t%i\t%i" % ("all",
                               all_ntotal,
                               len(all_genes),
                               len(all_categories)))

    return


def ReadGene2GOFromFile(infile, synonyms={}, obsolete={}):
    """reads GO mappings for all go_types from a
    file.

    If synonyms is given, goids in synynoms will be translated.
    Terms in *obsolete* will be discarded.

    returns two maps: gene2go maps genes to go categories
    and go2info maps go categories to information.
    """

    gene2gos = {}
    go2infos = {}
    c = E.Counter()

    for line in infile:
        if line[0] == "#":
            continue
        try:
            go_type, gene_id, goid, description, evidence = line[
                :-1].split("\t")
        except ValueError as msg:
            raise ValueError("parsing error in line '%s': %s" %
                             (line[:-1], msg))
        if go_type == "go_type":
            continue

        c.input += 1

        if goid in synonyms:
            c.synonyms += 1
            goid = synonyms[goid]

        if goid in obsolete:
            c.obsolete += 1
            continue

        gm = GOMatch(goid, go_type, description, evidence)
        gi = GOInfo(goid, go_type, description)
        if go_type not in gene2gos:
            gene2gos[go_type] = {}
            go2infos[go_type] = {}

        gene2go = gene2gos[go_type]
        go2info = go2infos[go_type]

        if gene_id not in gene2go:
            gene2go[gene_id] = []
        gene2go[gene_id].append(gm)
        go2info[goid] = gi
        c.output += 1

    E.debug("read gene2go assignments: %s" % str(c))

    return gene2gos, go2infos


def CountGO(gene2go):
    """count number of genes and go categories in mapping."""

    cats = collections.defaultdict(int)
    nmaps = 0
    for k, vv in list(gene2go.items()):
        for v in vv:
            nmaps += 1
            cats[v.mGOId] += 1

    return len(gene2go), len(cats), nmaps, cats


def removeCategories(gene2go, categories):
    '''remove all genes that map to *categories*.'''

    for k, vv in list(gene2go.items()):
        gene2go[k] = [v for v in vv if v.mGOId not in categories]


def countGOs(gene2gos):
    """return map of number of genes and go categories in mapping."""
    genes, goids = collections.defaultdict(int), collections.defaultdict(int)

    for cat, gene2go in gene2gos.items():
        for gene_id, vv in gene2go.items():
            genes[gene_id] += 1
            for v in vv:
                goids[v.mGOId] += 1
    return genes, goids


def ReadGeneLists(filename_genes, gene_pattern=None):
    """read gene lists from filename in matrix.

    returns a tuple (list of all genes, dictionary of gene lists)
    """

    if filename_genes == "-":
        infile = sys.stdin
    else:
        infile = IOTools.openFile(filename_genes, "r")

    headers, table = CSV.readTable(infile.readlines(), as_rows=False)

    if filename_genes != "-":
        infile.close()

    all_genes = table[0]

    # if there is only a single column, add a dummy column
    if len(table) == 1:
        table.append([1] * len(table[0]))
        headers.append("foreground")

    E.info("read %i genes from %s" % (len(all_genes), filename_genes))

    if gene_pattern:
        rx = re.compile(gene_pattern)
        all_genes = [rx.search(x).groups()[0] for x in all_genes]

    gene_lists = collections.OrderedDict()
    for header, col in zip(headers[1:], table[1:]):
        s = list(set([x for x, y in zip(all_genes, col) if y != "0"]))
        gene_lists[header] = set(s)

    return all_genes, gene_lists


def buildGO2Genes(gene2gos, ancestors=None):
    '''invert the dictionary genes2go.

    If ancestors is given, add missing ancestral information.
    '''

    go2genes = collections.defaultdict(set)
    for gene_id, terms in gene2gos.items():
        for term in terms:
            go2genes[term.mGOId].add(gene_id)
            if ancestors:
                for anc in ancestors[term.mGOId]:
                    go2genes[anc].add(gene_id)
    return go2genes


def GetCode(v):
    """return a code for over/underrepresentation."""

    if v.mRatio > 1.0:
        code = "+"
    elif v.mRatio < 1.0:
        code = "-"
    else:
        code = "?"
    return code


def convertGo2Goslim(options):
    """read gene list with GO assignments and convert to GO slim
    categories."""

    E.info("reading GO assignments from stdin")
    gene2gos, go2infos = ReadGene2GOFromFile(options.stdin)
    input_genes, input_goids = countGOs(gene2gos)

    #############################################################
    # read GO ontology from file
    assert options.filename_ontology, "please supply a GO ontology"
    E.info("reading ontology from %s" % (options.filename_ontology))

    infile = IOTools.openFile(options.filename_ontology)
    ontology = readOntology(infile)
    infile.close()

    go2infos = collections.defaultdict(dict)
    # substitute go2infos
    for go in list(ontology.values()):
        go2infos[go.mNameSpace][go.mId] = GOInfo(go.mId,
                                                 go_type=go.mNameSpace,
                                                 description=go.mName)

    E.info("reading GO assignments from %s" % options.filename_slims)
    go_slims = GetGOSlims(IOTools.openFile(options.filename_slims, "r"))

    if options.loglevel >= 1:
        v = set()
        for x in list(go_slims.values()):
            for xx in x:
                v.add(xx)
        E.info("read go slims from %s: go=%i, slim=%i" %
               (options.filename_slims,
                len(go_slims),
                len(v)))

    output_goids, output_genes = set(), set()
    noutput = 0
    options.stdout.write(
        "\t".join(("go_type", "gene_id", "go_id",
                   "description", "evidence")) + "\n")

    for category, gene2go in sorted(gene2gos.items()):
        gene2go = MapGO2Slims(gene2go, go_slims, ontology)
        for gene_id, values in sorted(gene2go.items()):
            output_genes.add(gene_id)
            for go in sorted(values, key=lambda x: x.mGOId):
                output_goids.add(go.mGOId)
                options.stdout.write("%s\t%s\t%s\t%s\t%s\n" %
                                     (go.mGOType,
                                      gene_id,
                                      go.mGOId,
                                      go.mDescription,
                                      "NA", ))
                noutput += 1

    E.info(
        ("ninput_genes=%i, ninput_goids=%i, noutput_gene=%i, "
         "noutput_goids=%i, noutput=%i") %
        (len(input_genes), len(input_goids),
         len(output_genes), len(output_goids),
         noutput))


def outputResults(outfile,
                  pairs,
                  go2info,
                  options,
                  fdrs=None,
                  samples=None,
                  gene2go=None,
                  foreground=None,
                  gene2name=None):
    '''output GO results to outfile.

    If foreground is given, output a list of gene identifiers in the
    foreground.

    If gene2name is given, output a columns with gene
    names (instead of identifiers)

    '''

    headers = ["code",
               "scount", "stotal", "spercent",
               "bcount", "btotal", "bpercent",
               "ratio",
               "pvalue", "pover", "punder",
               "goid", "category", "description"]

    if fdrs:
        headers += ["fdr"]

    if gene2go and foreground:
        headers += ['foreground']
        go2genes = buildGO2Genes(gene2go)
        if gene2name:
            headers += ['genes']

    if samples:
        headers += ["min", "max", "zscore", "mpover", "mpunder",
                    "nfdr_expected",
                    "CI95lower", "CI95upper"]

    outfile.write("\t".join(headers) + "\n")

    nselected = 0

    for k, v in sorted(pairs):

        code = GetCode(v)

        n = go2info.get(k, GOInfo())

        outfile.write("%s\t%s\t%s" % (code, str(v), n))

        if options.fdr:
            fdr = fdrs[k][0]
            outfile.write("\t%f" % fdr)

        if options.sample:

            if k in samples:
                s = samples[k]

            # calculate values for z-score
            if s.mStddev > 0:
                zscore = abs(
                    float(v.mSampleCountsCategory) - s.mMean) / s.mStddev
            else:
                zscore = 0.0

            outfile.write("\t%i\t%i\t%f\t%5.2e\t%5.2e\t%6.4f\t%6.4f\t%6.4f" %
                          (s.mMin,
                           s.mMax,
                           zscore,
                           min(s.mProbabilitiesOverRepresentation),
                           min(s.mProbabilitiesUnderRepresentation),
                           scipy.mean(s.mCounts),
                           scipy.stats.scoreatpercentile(s.mCounts, 5),
                           scipy.stats.scoreatpercentile(s.mCounts, 95),
                           ))

        if foreground:
            if k in go2genes:
                g = [x for x in go2genes[k] if x in foreground]
                if gene2name:
                    g = [gene2name.get(x, '?') for x in g]
                g = ";".join(sorted(g))
            else:
                g = ""
            outfile.write("\t%s" % g)

        outfile.write("\n")


def getSamples(gene2go, foreground, background, options, test_ontology,
               go2info):

    sample_size = options.sample
    # List of all minimum probabilities in simulation
    simulation_min_pvalues = []
    E.info("sampling: calculating %i samples: " % (sample_size))

    counts = {}
    prob_overs = {}
    prob_unders = {}

    samples = {}

    options.stdlog.write("# ")
    options.stdlog.flush()

    for x in range(sample_size):

        if options.loglevel >= 1:
            options.stdlog.write(".")
            options.stdlog.flush()

        # get shuffled array of genes from background
        sample_genes = random.sample(background, len(foreground))

        go_results = AnalyseGO(gene2go, sample_genes, background)

        pairs = list(go_results.mResults.items())

        for k, v in pairs:
            if k not in counts:
                counts[k] = []
                prob_overs[k] = []
                prob_unders[k] = []

            counts[k].append(v.mSampleCountsCategory)
            prob_overs[k].append(v.mProbabilityOverRepresentation)
            prob_unders[k].append(v.mProbabilityUnderRepresentation)

            simulation_min_pvalues.append(v.mPValue)

    if options.loglevel >= 1:
        sys.stdout.write("\n")
        sys.stdout.flush()

    E.info("sampling: sorting %i P-Values" % len(simulation_min_pvalues))

    simulation_min_pvalues.sort()
    simulation_min_pvalues = numpy.array(simulation_min_pvalues)

    samples = {}

    if options.output_filename_pattern:
        filename = options.output_filename_pattern % {
            'go': test_ontology, 'section': "samples"}
        E.info("sampling results go to %s" % filename)
        outfile = IOTools.openFile(filename, "w", create_dir=True)
    else:
        outfile = sys.stdout

    outfile.write("\t".join(("goid", "min", "max", "mean", "median", "stddev",
                             "CI95lower", "CI95upper",
                             "pover", "punder", "goid",
                             "category", "description")) + "\n")
    for k in sorted(list(counts.keys())):

        c = counts[k]

        prob_overs[k].sort()
        prob_unders[k].sort()

        s = GOSample(min(c),
                     max(c),
                     scipy.mean(c),
                     numpy.std(c),
                     numpy.array(prob_overs[k]),
                     numpy.array(prob_unders[k]),
                     counts[k])

        samples[k] = s

        outfile.write("%s\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" %
                      (k,
                       min(c),
                       max(c),
                       scipy.mean(c),
                       scipy.median(c),
                       numpy.std(c),
                       scipy.stats.scoreatpercentile(c, 5),
                       scipy.stats.scoreatpercentile(c, 95),
                       min(prob_overs[k]),
                       min(prob_unders[k]),
                       go2info[k]))

    if options.output_filename_pattern:
        outfile.close()

    return samples, simulation_min_pvalues


def computeFDRs(go_results,
                foreground,
                background,
                options,
                test_ontology,
                gene2go,
                go2info):

    pairs = sorted(go_results.mResults.items())

    E.info("calculating the FDRs using method `%s`" % options.qvalue_method)

    samples = None

    observed_min_pvalues = [min(x[1].mProbabilityOverRepresentation,
                                x[1].mProbabilityUnderRepresentation) for x in pairs]

    fdrs = {}

    method = options.qvalue_method

    if options.qvalue_method == "storey":

        # compute fdr via Storey's method
        try:
            fdr_data = Stats.doFDR(observed_min_pvalues)

        except ValueError as msg:
            E.warn("failure in q-value computation: %s" % msg)
            E.warn("reverting to Bonferroni correction")
            method = "bonf"
            fdr_data = Stats.FDRResult()
            l = float(len(observed_min_pvalues))
            fdr_data.mQValues = [min(1.0, x * l) for x in observed_min_pvalues]

        for pair, qvalue in zip(pairs, fdr_data.mQValues):
            fdrs[pair[0]] = (qvalue, 1.0, 1.0)

    elif options.qvalue_method == "empirical":
        assert options.sample > 0, "requiring a sample size of > 0"

        #######################################################################
        # sampling
        # for each GO-category:
        # get maximum and minimum counts in x samples -> calculate minimum/maximum significance
        # get average and stdev counts in x samples -> calculate z-scores for
        # test set
        samples, simulation_min_pvalues = getSamples(gene2go,
                                                     foreground,
                                                     background,
                                                     options,
                                                     test_ontology,
                                                     go2info)

        # compute P-values from sampling
        observed_min_pvalues.sort()
        observed_min_pvalues = numpy.array(observed_min_pvalues)

        sample_size = options.sample

        for k, v in pairs:

            if k in samples:
                s = samples[k]
            else:
                raise KeyError("category %s not in samples" % k)

            # calculate values for z-score
            if s.mStddev > 0:
                zscore = abs(
                    float(v.mSampleCountsCategory) - s.mMean) / s.mStddev
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
            pvalue = v.mPValue

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
                fdr = 1.0

            fdrs[k] = (fdr, a, b)
    else:
        qvalues = R['p.adjust'](
            observed_min_pvalues, method=options.qvalue_method)
        fdr_data = Stats.FDRResult()
        fdr_data.mQValues = list(qvalues)
        for pair, qvalue in zip(pairs, fdr_data.mQValues):
            fdrs[pair[0]] = (qvalue, 1.0, 1.0)

    return fdrs, samples, method


def getFileName(options, **kwargs):
    '''return a filename

    Placeholders in filename are string-substituted with the
    dictionary in kwargs.
    '''
    if options.output_filename_pattern:
        filename = options.output_filename_pattern % kwargs
        E.info("output for section '%s' go to %s" %
               (kwargs.get("section", "unknown"), filename))
        outfile = IOTools.openFile(filename, "w", create_dir=True)
    else:
        outfile = options.stdout

    return outfile


def buildMatrix(results, valuef, dtype=numpy.float, default=0):
    '''build a matrix from a field in *results*

    The value stored in the matrix is accessed via *valuef*.
    '''

    row_headers = [set([x[0] for x in y]) for y in results]
    row_headers = sorted(list(row_headers[0].union(*row_headers[1:])))
    map_row = dict(list(zip(row_headers, list(range(len(row_headers))))))
    matrix = numpy.zeros((len(row_headers), len(results)), dtype=dtype)
    if default != 0:
        matrix[:] = default

    for col, pairs in enumerate(results):
        for row, v in pairs:
            try:
                matrix[map_row[row]][col] = valuef(v)
            except ValueError:
                # ignore errors for log(0)
                pass

    return matrix, row_headers


def selectSignificantResults(pairs, fdrs, options):
    '''select a set of significant results.
    '''

    filtered_pairs = []
    for k, v in pairs:

        is_ok = False

        pvalue = v.mPValue

        if options.fdr:
            (fdr, expfpos, pos) = fdrs[k]
            if fdr < options.threshold:
                is_ok = True
        else:
            if pvalue < options.threshold:
                is_ok = True

        if is_ok:
            filtered_pairs.append((k, v))

    return filtered_pairs


def outputMultipleGeneListResults(results,
                                  all_genelists_with_results,
                                  test_ontology,
                                  go2info,
                                  options,
                                  section):
    '''select a set of significant results.
    '''

    col_headers = all_genelists_with_results

    if len(results) == 0:
        E.warn('no significant results - no matrices output')
        return

    assert len(col_headers) == len(results)

    def _output(section, subsection, valuef, dtype):

        # fold change matrix
        matrix, row_headers = buildMatrix(results,
                                          valuef=valuef,
                                          dtype=dtype)

        outfile = getFileName(options,
                              go=test_ontology,
                              section=section,
                              set='%s_all' % subsection)

        IOTools.writeMatrix(
            outfile, matrix, row_headers, col_headers, row_header="category")

        outfile = getFileName(options,
                              go=test_ontology,
                              section=section,
                              set='%s_alldesc' % subsection)

        IOTools.writeMatrix(outfile, matrix,
                            ["%s:%s" % (x, go2info[x].mDescription)
                             for x in row_headers],
                            col_headers, row_header="category")

    _output('l2fold', section,
            valuef=lambda x: math.log(x.mRatio + 0.00000001, 2),
            dtype=numpy.float)

    _output('l10pvalue', section,
            valuef=lambda x: int(-10 * math.log(x.mPValue, 10)),
            dtype=numpy.int)

    _output('l10qvalue', section,
            valuef=lambda x: int(-10 * math.log(x.mQValue, 10)),
            dtype=numpy.int)


def pairwiseGOEnrichment(results_per_genelist, labels, test_ontology, go2info,
                         options):
    '''compute pairwise enrichment between sets.

    The purpose of this method is to find if there are categories that are differently enriched
    in a pair of gene lists.

    The appropriate test here is the Chi-Squared test.

    The assumption is that the background set is the same in all gene lists.

    The workflow is thus::

       for each combination of two gene lists:
           for each GO category:
               get counts in foreground, total counts of foreground
               compute chi-square enrichment output
               save P-value
           apply fdr - output significant differences.
    '''

    dicts = [dict(x) for x in results_per_genelist]

    PairResult = collections.namedtuple("PairResult",
                                        "goid set1 set2 counts1 total1 pvalue1 qvalue1 counts2 total2 pvalue2 qvalue2 pvalue qvalue description")

    outfile = getFileName(options,
                          go=test_ontology,
                          section='summary',
                          set="pairs")

    outfile.write(
        "set1\tset2\ttotal1\ttotal2\tshared\tskipped\ttested\tsignificant\tinsignificant\n")

    results = []

    total = len(dicts) * (len(dicts) - 1) / 2

    iteration = 0

    min_observed_counts = options.pairs_min_observed_counts

    for x, genelist1 in enumerate(sorted(dicts)):

        x_go_categories = set(genelist1.keys())
        for y, genelist2 in enumerate(sorted(dicts[:x])):

            iteration += 1
            if iteration % 10 == 0:
                E.info("iteration: %i/%i (%5.2f%%)" %
                       (iteration, total, 100.0 * iteration / total))

            y_go_categories = set(genelist2.keys())

            shared = x_go_categories.intersection(y_go_categories)

            c = E.Counter()

            for category in shared:
                c.shared += 1
                xx = genelist1[category]
                yy = genelist2[category]

                # discard all tests with few observations in the observed
                # counts
                if xx.mSampleCountsCategory < min_observed_counts and yy.mSampleCountsCategory < min_observed_counts:
                    c.skipped += 1
                    continue

                observed = (xx.mSampleCountsCategory, yy.mSampleCountsCategory)

                aa, bb, cc, dd = \
                    (xx.mSampleCountsCategory,
                     yy.mSampleCountsCategory,
                     xx.mSampleCountsTotal - xx.mSampleCountsCategory,
                     yy.mSampleCountsTotal - yy.mSampleCountsCategory)

                if cc == dd == 0:
                    c.skipped += 1
                    continue

                c.tested += 1

                fisher, pvalue = scipy.stats.fisher_exact(numpy.array(
                    ((aa, bb),
                     (cc, dd))))

                if pvalue < 0.05:
                    c.significant_pvalue += 1
                else:
                    c.insignificant_pvalue += 1

                results.append(PairResult._make((category,
                                                 labels[x],
                                                 labels[y],
                                                 xx.mSampleCountsCategory,
                                                 xx.mSampleCountsTotal,
                                                 xx.mPValue,
                                                 xx.mQValue,
                                                 yy.mSampleCountsCategory,
                                                 yy.mSampleCountsTotal,
                                                 yy.mPValue,
                                                 yy.mQValue,
                                                 pvalue,
                                                 1.0,
                                                 go2info[category].mDescription)))

            outfile.write("\t".join(map(str,
                                        (labels[x], labels[y],
                                         len(x_go_categories),
                                         len(y_go_categories),
                                         c.shared,
                                         c.skipped,
                                         c.tested,
                                         c.significant_pvalue,
                                         c.insignicant_pvalue))) + "\n")
    if options.output_filename_pattern:
        outfile.close()

    if options.fdr:
        pvalues = [x.pvalue for x in results]

        if options.qvalue_method == "storey":

            # compute fdr via Storey's method
            try:
                fdr_data = Stats.doFDR(pvalues)

            except ValueError as msg:
                E.warn("failure in q-value computation: %s" % msg)
                E.warn("reverting to Bonferroni correction")
                method = "bonf"
                fdr_data = Stats.FDRResult()
                l = float(len(pvalues))
                fdr_data.mQValues = [min(1.0, x * l) for x in pvalues]

            qvalues = fdr_data.mQValues
        else:
            qvalues = R['p.adjust'](pvalues, method=options.qvalue_method)

        # update qvalues
        results = [x._replace(qvalue=y) for x, y in zip(results, qvalues)]

    outfile = getFileName(options,
                          go=test_ontology,
                          section='pairs',
                          set="pairs")

    outfile.write("\t".join(PairResult._fields) + "\n")
    for result in results:
        outfile.write("\t".join(map(str, result)) + "\n")

    if options.output_filename_pattern:
        outfile.close()
