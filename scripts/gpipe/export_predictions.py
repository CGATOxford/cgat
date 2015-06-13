##########################################################################
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
##########################################################################
'''
gpipe/export_predictions.py - 
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

   python gpipe/export_predictions.py --help

Type::

   python gpipe/export_predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import time
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pgdb
import CGAT.IndexedFasta as IndexedFasta

USAGE = """python %s [OPTIONS] < in > out

output gene predictions for AAA submission.

Version: $Id: gpipe/export_predictions.py 2781 2009-09-10 11:33:14Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]


##########################################################################
##########################################################################
##########################################################################
# Class to select entries for tracks (genes/cds/rnas)
##########################################################################


class Selector:

    def __init__(self, dbhandle, options):

        self.mOptions = options

        self.mDbHandle = dbhandle

        # pre-processing: set table names
        self.mParameters = {}
        self.mParameters['tablename_predictions'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_predictions
        self.mParameters['tablename_genes'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_genes
        self.mParameters['tablename_redundant'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_redundant
        self.mParameters['tablename_contigs'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_contigs
        self.mParameters['tablename_exons'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_exons
        self.mParameters['tablename_quality'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_quality
        self.mParameters['tablename_disruptions'] = self.mOptions.schema + \
            "." + self.mOptions.tablename_disruptions

    def readFilter(self):
        """read a list of transcripts/genes to be filtered."""
        if self.mOptions.filename_filter:
            self.mFilterList, errors = IOTools.ReadList(
                open(self.mOptions.filename_filter), map_function=int)

        if len(self.mFilterList) == 0:
            raise ValueError("filter is empty, no predictions will be output.")

    def addFilter(self):
        """add options to select statement (in the where part). Start with an 'AND'."""
        return ""

    def addTables(self):
        """add tables to the tables for (for filtering). Start with a ','."""
        return ""

    def getRows(self):
        """execute a select statement. returns a list of hashs."""

        statement = self.getStatement()

        cc = self.mDbHandle.cursor()
        cc.execute(statement)

        descriptions = cc.description

        rows = []
        fields = range(len(descriptions))
        while 1:
            r = cc.fetchone()
            if not r:
                break
            rr = {}
            for f in fields:
                rr[descriptions[f][0]] = r[f]
            rows.append(rr)
        cc.close()

        return rows

# ------------------------------------------------------------------------


class SelectorDisruptions(Selector):

    def __init__(self, dbhandle, options):

        Selector.__init__(self, dbhandle, options)

    def getStatement(self, type):

        statement = """
        SELECT
        p.sbjct_token AS sbjct_token,
        p.sbjct_strand AS sbjct_strand,
        g.gene_id AS gene_id,
        p.prediction_id as prediction_id, 
        MIN(CASE WHEN p.sbjct_strand = '+' THEN d.genome_from+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-d.genome_to+c.start END)+1 AS feature_from, 
        MAX(CASE WHEN p.sbjct_strand = '+' THEN d.genome_to+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-d.genome_from+c.start END) AS feature_to,
        cds_from AS cds_from,
        cds_to AS cds_to
        FROM
        %(tablename_predictions)s AS p,
        %(tablename_genes)s AS g,
        %(tablename_contigs)s AS c,
        %(tablename_disruptions)s AS d
        WHERE \
        g.prediction_id = p.prediction_id AND
        c.sbjct_token = p.sbjct_token AND
        d.prediction_id = p.prediction_id AND
        d.type = '%%s'
        """ % (self.mParameters)

        statement = statement % (type)

        statement += """
        GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id,
        d.cds_from, d.cds_to"""

        statement += " ORDER BY p.sbjct_token, feature_from"

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "statement for selecting disruptions:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

        return statement

    def getRows(self, type):
        """execute a select statement. returns a list of hashs."""

        statement = self.getStatement(type)

        cc = self.mDbHandle.cursor()
        cc.execute(statement)

        descriptions = cc.description

        rows = []
        fields = range(len(descriptions))
        while 1:
            r = cc.fetchone()
            if not r:
                break
            rr = {}
            for f in fields:
                rr[descriptions[f][0]] = r[f]
            rows.append(rr)
        cc.close()

        return rows


# ------------------------------------------------------------------------
class SelectorMRNA(Selector):

    def __init__(self, dbhandle, options):

        Selector.__init__(self, dbhandle, options)

    def getStatement(self):

        statement = """
        SELECT
        p.sbjct_token AS sbjct_token,
        p.sbjct_strand AS sbjct_strand,
        CASE WHEN g.gene_id = '0' THEN g.gene_id || '_' || p.prediction_id ELSE CAST( g.gene_id AS TEXT) END AS gene_id,
        p.prediction_id AS prediction_id,
        MIN(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS feature_from, 
        MAX(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END) AS feature_to,
        p.query_token AS query_token,
        p.query_coverage AS query_coverage,
        p.pidentity AS percent_identity,
        p.nframeshifts AS nframeshifts,        
        p.nstopcodons AS nstopcodons,
        q.class AS class
        FROM
        %(tablename_predictions)s AS p,
        %(tablename_genes)s AS g,
        %(tablename_exons)s AS e,
        %(tablename_quality)s AS q,
        %(tablename_contigs)s AS c
        %%s
        WHERE g.prediction_id = p.prediction_id AND
        e.prediction_id = p.prediction_id AND
        q.prediction_id = p.prediction_id AND
        c.sbjct_token = p.sbjct_token AND
        e.genome_exon_from > 0 AND e.genome_exon_to > 0
        %%s
        """ % (self.mParameters)

        statement = statement % (self.addTables(), self.addFilter())

        statement += """
        GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id,
        p.query_token, p.query_coverage, p.pidentity, p.nframeshifts, p.nstopcodons, q.class"""

        statement += " ORDER BY p.sbjct_token, feature_from"

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "statement for selecting mRNAs:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

        return statement

# ------------------------------------------------------------------------


class SelectorGenes(Selector):

    def __init__(self, dbhandle, options):

        Selector.__init__(self, dbhandle, options)

    def getStatement(self):
        """get statement.

        Gene ids for gene 0 are modified such that the gene_id will be
        0_prediction_id.
        """

        statement = """
        SELECT
        p.sbjct_token AS sbjct_token,
        p.sbjct_strand AS sbjct_strand,
        CASE WHEN g.gene_id = '0' THEN g.gene_id || '_' || p.prediction_id ELSE CAST (g.gene_id AS TEXT) END AS gene_id,
        MIN(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS feature_from, 
        MAX(CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END) AS feature_to
        FROM
        %(tablename_predictions)s AS p,
        %(tablename_genes)s AS g,
        %(tablename_exons)s AS e,
        %(tablename_contigs)s AS c
        %%s
        WHERE g.prediction_id = p.prediction_id AND
        e.prediction_id = p.prediction_id AND
        c.sbjct_token = p.sbjct_token AND
        e.genome_exon_from > 0 AND e.genome_exon_to > 0
        %%s
        """ % (self.mParameters)

        statement = statement % (self.addTables(), self.addFilter())

        statement += " GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id"
        statement += " ORDER BY p.sbjct_token, feature_from"

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "statement for selecting genes:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

        return statement

# ------------------------------------------------------------------------


class SelectorCDS(Selector):

    def __init__(self, dbhandle, options):

        Selector.__init__(self, dbhandle, options)

    def getStatement(self):
        """there might be several entries for one exons if an exon is joined.
        """

        statement = """
        SELECT DISTINCT
        p.sbjct_token AS sbjct_token,
        p.sbjct_strand AS sbjct_strand,
        CASE WHEN g.gene_id = '0' THEN g.gene_id || '_' || p.prediction_id ELSE CAST( g.gene_id AS TEXT) END AS gene_id,
        p.prediction_id AS prediction_id,
        (CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_from+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_to+c.start END)+1 AS feature_from, 
        CASE WHEN p.sbjct_strand = '+' THEN e.genome_exon_to+c.start 
        WHEN p.sbjct_strand = '-' THEN c.size-e.genome_exon_from+c.start END AS feature_to,
        BOOL_OR(e.is_ok) AS is_ok,
        MAX(e.pidentity) AS percent_identity,
        e.exon_id AS exon_id,
        e.exon_frame AS frame
        FROM
        %(tablename_predictions)s AS p,
        %(tablename_genes)s AS g,
        %(tablename_exons)s AS e,
        %(tablename_contigs)s AS c
        %%s
        WHERE g.prediction_id = p.prediction_id AND
        e.prediction_id = p.prediction_id AND
        c.sbjct_token = p.sbjct_token AND
        e.genome_exon_from > 0 AND e.genome_exon_to > 0
        %%s
        """ % (self.mParameters)

        statement = statement % (self.addTables(), self.addFilter())

        statement += " GROUP BY p.sbjct_token, p.sbjct_strand, g.gene_id, p.prediction_id, feature_from, feature_to, e.exon_id, e.exon_frame"
        statement += " ORDER BY p.sbjct_token, feature_from"

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "statement for selecting cds:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

        return statement

# ------------------------------------------------------------------------


class SelectorContigs(Selector):

    """select contigs."""

    def __init__(self, dbhandle, options):
        Selector.__init__(self, dbhandle, options)

    def getStatement(self):
        return """
        SELECT c.sbjct_token as sbjct_token, 1 as feature_from, c.size AS feature_to
        FROM
        %(tablename_contigs)s AS c    
        """ % (self.mParameters)

##########################################################################
##########################################################################
##########################################################################
# Class to select predictions (genes/predictions)
##########################################################################


class SelectorPredictions(Selector):

    def __init__(self, dbhandle, options):
        Selector.__init__(self, dbhandle, options)

    def getStatement(self):
        """get statement for selecting predictions.
        """

        statement = """
        SELECT DISTINCT
        p.prediction_id AS prediction_id,
        CASE WHEN g.gene_id = '0' THEN g.gene_id || '_' || p.prediction_id ELSE CAST (g.gene_id AS TEXT) END AS gene_id
        FROM
        %(tablename_predictions)s AS p,
        %(tablename_genes)s AS g
        %%s
        WHERE g.prediction_id = p.prediction_id 
        %%s
        """ % (self.mParameters)

        statement = statement % (self.addTables(), self.addFilter())

        statement += " GROUP BY g.gene_id, p.prediction_id"

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "statement for selecting:\n%s\n%s\n%s\n" % ("#" * 20, statement, "#" * 20))

        return statement

    def getRows(self):

        statement = self.getStatement()

        cc = self.mDbHandle.cursor()
        cc.execute(statement)

        descriptions = cc.description

        predictions, genes = {}, {}

        for prediction_id, gene_id in cc.fetchall():
            predictions[str(prediction_id)] = str(gene_id)
            genes[str(gene_id)] = True

        if self.mOptions.loglevel >= 2:
            self.mOptions.stdlog.write(
                "# retrieved %i predictions in %i genes.\n" % (len(predictions), len(genes)))

        return predictions, genes

##########################################################################
# Subselect classes
# All predictions


class SelectorPredictionsRaw(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

##########################################################################
# Filtered predictions: all non-redundant predictions


class SelectorPredictionsFiltered(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s.%s AS r" % (options.schema, options.tablename_redundant)

    def addFilter(self):
        return """AND g.gene_id > 0
    AND p.prediction_id = r.rep_prediction_id"""

##########################################################################
# Clean predictions


class SelectorPredictionsClean(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s.%s AS r,%s.%s AS q" % (options.schema, options.tablename_redundant,
                                           options.schema, options.tablename_quality)

    def addFilter(self):
        return """AND g.gene_id > 0 
    AND p.prediction_id = r.rep_prediction_id
    AND p.prediction_id = q.prediction_id
    AND q.class IN ('%s')
    """ % "','".join(self.mOptions.quality_clean_set)

##########################################################################
# Clean predictions


class SelectorPredictionsOptic(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s.%s AS r,%s.%s AS q" % (options.schema, options.tablename_redundant,
                                           options.schema, options.tablename_quality)

    def addFilter(self):
        return """AND g.gene_id > 0 
    AND p.prediction_id = r.rep_prediction_id
    AND p.prediction_id = q.prediction_id
    AND q.class IN ('%s')
    """ % "','".join(self.mOptions.quality_optic_set)

##########################################################################
# Predictions with orthologs


class SelectorPredictionsOrthologs(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s AS oo" % options.tablename_orthologs

    def addFilter(self):
        return "AND g.gene_id > 0 AND oo.prediction_id1 = p.prediction_id AND oo.schema1='%s' " % (self.mOptions.schema)

##########################################################################
# Predictions with a minimum ds


class SelectorPredictionsMaxKs(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s.%s AS k" % (options.schema, self.mOptions.tablename_kaks)

    def addFilter(self):
        return "AND g.gene_id > 0 AND k.prediction_id = p.prediction_id AND k.ds <= %f " % (self.mOptions.max_ks)

##########################################################################
# Predictions with a minimum ds


class SelectorPredictionsMaxKsQuality(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)

    def addTables(self):
        return ",%s.%s AS k, %s.%s AS q" % (options.schema, self.mOptions.tablename_kaks,
                                            options.schema, self.mOptions.tablename_quality)

    def addFilter(self):
        return """AND g.gene_id > 0
    AND k.prediction_id = p.prediction_id
    AND p.prediction_id = q.prediction_id 
    AND ( (q.class IN ('CG')) OR
    (q.class IN ('PG', 'SG', 'CF', 'PF') AND k.ds <= %f))
    """ % (self.mOptions.max_ks)

##########################################################################
# Select predictions from a list of predictions


class SelectorPredictionsFilterPredictions(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)
        self.readFilter()

    def addFilter(self):
        return """AND p.prediction_id IN ('%s')""" % ("','".join(map(str, self.mFilterList)))


class SelectorPredictionsFilterGenes(SelectorPredictions):

    def __init__(self, dbhandle, options):
        SelectorPredictions.__init__(self, dbhandle, options)
        self.readFilter()

    def addFilter(self):
        return """AND g.gene_id IN ('%s')""" % ("','".join(map(str, self.mFilterList)))


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
class Exporter:

    def __init__(self,
                 selector_predictions,
                 selector_genes, selector_mrnas, selector_cds,
                 options, do_map=False,
                 selector_contigs=None,
                 selector_disruptions=None):

        self.mDoMap = do_map

        self.mOptions = options

        self.mSelectorGenes = selector_genes
        self.mSelectorMRNAs = selector_mrnas
        self.mSelectorCDS = selector_cds
        self.mSelectorContigs = selector_contigs
        self.mSelectorPredictions = selector_predictions
        self.mSelectorDisruptions = selector_disruptions

        self.mMethod = "gpipe"

        self.mPredictions, self.mGenes = selector_predictions.getRows()

        # filtering is done via the selectors,
        # this silences outputs for certain tracks.
        self.mWriteMRNATrack = True
        self.mWriteGeneTrack = True
        self.mWriteCDSTrack = True
        self.mWriteExonTrack = False
        self.mWriteStopCodon = False
        self.mWriteStartCodon = False

        self.mExtraTracks = None

    # ------------------------------------------------------------------------
    def addTrack(self, track):
        """add a track to extra tracks."""
        self.mExtraTracks.append(track)

    # ------------------------------------------------------------------------
    def setMethod(self, method):
        self.mMethod = method

    # ------------------------------------------------------------------------
    def getGeneIdentifier(self, gene_id):
        return str(gene_id)

    def getMRNAIdentifier(self, prediction_id, gene_id):
        return str(prediction_id)

    # ------------------------------------------------------------------------
    def getInfoGene(self, gene_id, row):
        """get gene info."""
        return "Id=%s" % gene_id

    # ------------------------------------------------------------------------
    def getInfoMRNA(self, prediction_id, gene_id, row):
        """get mrna info."""
        return ";".join(("Id=%s" % prediction_id,
                         "Parent=%s" % gene_id))

    # ------------------------------------------------------------------------
    def getInfoCDS(self, prediction_id, gene_id, row):
        """get cds info."""
        return ";".join(("Id=%s" % prediction_id,
                         "Parent=%s" % gene_id))

    # ------------------------------------------------------------------------
    def getInfoDisruption(self, prediction_id, gene_id, row):
        """get info for stop codons."""
        return ";".join(("Id=%s" % prediction_id,
                         "cds_from=%i" % (row['cds_from']),
                         "cds_to=%i" % (row['cds_to'])))

    # ------------------------------------------------------------------------
    def writeGeneTrack(self,
                       outfile,
                       current_id):

        rr = self.mSelectorGenes.getRows()

        map_gene_id = {}

        for row in rr:

            gene_id = str(row['gene_id'])

            if gene_id not in self.mGenes:
                continue

            if not self.mDoMap:
                current_id = gene_id
            else:
                current_id += 1

            if gene_id not in map_gene_id:
                map_gene_id[gene_id] = self.getGeneIdentifier(current_id)

            if self.mWriteGeneTrack:
                outfile.write("\t".join((row['sbjct_token'],
                                         self.mMethod,
                                         "gene",
                                         str(row['feature_from']),
                                         str(row['feature_to']),
                                         ".",
                                         row['sbjct_strand'],
                                         ".",
                                         self.getInfoGene(map_gene_id[gene_id], row))) + "\n")

        return current_id, map_gene_id

    # ------------------------------------------------------------------------
    def writeMRNATrack(self,
                       outfile,
                       current_id, map_gene_id):

        rr = self.mSelectorMRNAs.getRows()

        map_prediction_id = {}
        for row in rr:

            gene_id = str(row['gene_id'])
            prediction_id = str(row['prediction_id'])

            if prediction_id not in self.mPredictions:
                continue

            if gene_id not in self.mGenes:
                continue

            if gene_id not in map_gene_id:
                continue

            if not self.mDoMap:
                current_id = prediction_id
            else:
                current_id += 1

            if prediction_id not in map_prediction_id:
                map_prediction_id[prediction_id] = self.getMRNAIdentifier(
                    current_id, map_gene_id[gene_id])

            if self.mWriteMRNATrack:
                outfile.write("\t".join((row['sbjct_token'],
                                         self.mMethod,
                                         "mRNA",
                                         str(row['feature_from']),
                                         str(row['feature_to']),
                                         ".",
                                         row['sbjct_strand'],
                                         ".",
                                         self.getInfoMRNA(map_prediction_id[prediction_id], map_gene_id[gene_id], row))) + "\n")

        return current_id, map_prediction_id

    # ------------------------------------------------------------------------
    def writeDisruptions(self,
                         outfile,
                         map_gene_id, map_prediction_id,
                         type="stop", track="stop"):

        rr = self.mSelectorDisruptions.getRows(type)

        for row in rr:

            gene_id = str(row['gene_id'])
            prediction_id = str(row['prediction_id'])

            if prediction_id not in self.mPredictions:
                continue

            if gene_id not in self.mGenes:
                continue

            if (gene_id not in map_gene_id) or \
               (prediction_id not in map_prediction_id):
                continue

            outfile.write("\t".join((row['sbjct_token'],
                                     self.mMethod,
                                     track,
                                     str(row['feature_from']),
                                     str(row['feature_to']),
                                     ".",
                                     row['sbjct_strand'],
                                     ".",
                                     self.getInfoDisruption(map_prediction_id[prediction_id],
                                                            map_gene_id[
                                                                gene_id],
                                                            row))) + "\n")

    # ------------------------------------------------------------------------
    def writeCDSTrack(self,
                      outfile,
                      current_id,
                      map_gene_id, map_prediction_id):

        rr = self.mSelectorCDS.getRows()

        for row in rr:

            gene_id = str(row['gene_id'])
            prediction_id = str(row['prediction_id'])

            if prediction_id not in self.mPredictions:
                continue

            if gene_id not in self.mGenes:
                continue

            if gene_id not in map_gene_id:
                continue

            if prediction_id not in map_prediction_id:
                continue

            if self.mWriteCDSTrack:
                outfile.write("\t".join((row['sbjct_token'],
                                         self.mMethod,
                                         "CDS",
                                         str(row['feature_from']),
                                         str(row['feature_to']),
                                         ".",
                                         row['sbjct_strand'],
                                         str(row['frame']),
                                         self.getInfoCDS(map_prediction_id[prediction_id],
                                                         map_gene_id[
                                                             gene_id],
                                                         row))) + "\n")

                # if desired, write an exon for each CDS
                if self.mWriteExonTrack:
                    outfile.write("\t".join((row['sbjct_token'],
                                             self.mMethod,
                                             "exon",
                                             str(row['feature_from']),
                                             str(row['feature_to']),
                                             ".",
                                             row['sbjct_strand'],
                                             ".",
                                             self.getInfoCDS(map_prediction_id[prediction_id],
                                                             map_gene_id[
                                                                 gene_id],
                                                             row))) + "\n")

    # ------------------------------------------------------------------------
    def writeSequenceRegions(self,
                             outfile, dbhandle):

        statement = """
        SELECT c.sbjct_token, 1, c.size
        FROM
        %(tablename_contigs)s AS c    
        """ % (self.mParameters)

        cc = dbhandle.cursor()
        cc.execute(statement)
        rr = cc.fetchall()
        cc.close()

        for row in rr:
            outfile.write(" ".join(("##sequence-region",
                                    row['sbjct_token'],
                                    str(row['feature_from']),
                                    str(row['feature_to']),
                                    )) + "\n")

    # ------------------------------------------------------------------------
    def writeContigs(self,
                     outfile):

        if not self.mSelectorContigs:
            return

        rr = self.mSelectorContigs.getRows()

        for row in rr:
            outfile.write("\t".join((row['sbjct_token'],
                                     self.mMethod,
                                     "Contig",
                                     str(row['feature_from']),
                                     str(row['feature_to']),
                                     ".",
                                     ".",
                                     ".",
                                     "Contig %s" % row['sbjct_token']
                                     )) + "\n")

    # ------------------------------------------------------------------------
    def writeHeader(self, outfile):
        """write Header (default = no header)."""
        return

    # ------------------------------------------------------------------------
    def writeGFF(self, outfile):

        # start numbering from the first id
        current_id = self.mOptions.first_id

        # output header
        self.writeHeader(outfile)

        # write sequence regions
        if "regions" in self.mOptions.tracks:
            self.writeSequenceRegions(outfile)

        # write tracks
        if "contigs" in self.mOptions.tracks:
            self.writeContigs(outfile)

        map_gene_id, map_prediction_id = None, None

        # write the various tracks major tracks. These are necessary to map
        # predictions to genes and renumber them. Output might be silenced.
        if "genes" in self.mOptions.tracks:
            current_id, map_gene_id = self.writeGeneTrack(outfile,
                                                          current_id)

        if self.mOptions.loglevel > 3:
            print "# writing the following %i genes: " % len(map_gene_id)
            print "#", map_gene_id

        if "mrnas" in self.mOptions.tracks:
            current_id, map_prediction_id = self.writeMRNATrack(outfile,
                                                                current_id,
                                                                map_gene_id)

        if "cds" in self.mOptions.tracks:
            self.writeCDSTrack(outfile, current_id,
                               map_gene_id, map_prediction_id)

        if map_gene_id and map_prediction_id:
            self.mOptions.stdlog.write("# genome=%s, genes=%i, transcripts=%i\n" %
                                       (self.mOptions.schema, len(map_gene_id), len(map_prediction_id)))
            self.mOptions.stdlog.flush()

        if "stops" in self.mOptions.tracks:
            self.writeDisruptions(outfile, map_gene_id, map_prediction_id,
                                  type="stop", track="stop")
            self.writeDisruptions(outfile, map_gene_id, map_prediction_id,
                                  type="split-stop", track="stop")

        if "frameshifts" in self.mOptions.tracks:
            self.writeDisruptions(outfile, map_gene_id, map_prediction_id,
                                  type="frameshift", track="frameshift")

        return map_gene_id, map_prediction_id

    # ------------------------------------------------------------------------
    def writeMap(self, outfile, category, map_id):
        """write map between genes/predictions to new ids."""

        outfile = open(self.mOptions.filename_pattern_map %
                       (self.mOptions.assembly_id, category), "w")
        ids = map_id.keys()
        ids.sort()
        for id in ids:
            outfile.write("%i\t%s\n" % (id, map_id[id]))
        outfile.close()

##########################################################################
##########################################################################
##########################################################################
##########################################################################
# Exporter subclasses
##########################################################################

# ------------------------------------------------------------------------


class ExporterAAA(Exporter):

    def __init__(self,
                 selector_predictions,
                 selector_genes, selector_mrnas, selector_cds,
                 options,
                 selector_contigs=None,
                 selector_disruptions=None):

        Exporter.__init__(self,
                          selector_predictions,
                          selector_genes, selector_mrnas, selector_cds,
                          options, do_map=True,
                          selector_contigs=selector_contigs,
                          selector_disruptions=selector_disruptions)

        self.mMethod = self.mOptions.group_id + "_" + self.mOptions.method_id

    def getGeneIdentifier(self, gene_id):
        return "%s_%s_%s_%s" % (self.mOptions.species_id, self.mOptions.group_id, self.mOptions.method_id, gene_id)

    def getMRNAIdentifier(self, prediction_id, gene_id):
        return "%s_%s_%s_%s" % (self.mOptions.species_id, self.mOptions.group_id, self.mOptions.method_id, prediction_id)

    def writeHeader(self, outfile):

        outfile.write("""##gff-version   3
#species: %(species)s
#assembly-id: %(assembly)s
#annotation-group-id: %(group_id)s_%(method_id)s
#algorithm: %(group_id)s_%(method_id)s = Oxford gene prediction pipeline based on the program exonerate.
#templates: Gene models from FlyBase 4.2.1 obtained from ENSEMBL version 37
#authors: Andreas Heger and Chris Ponting at firstname.lastname@anat.ox.ac.uk
#date: %(date)s
""" % { 'species': self.mOptions.species,
            'assembly': self.mOptions.assembly_id,
            'group_id': self.mOptions.group_id,
            'method_id': self.mOptions.method_id,
            'date': time.strftime("%Y%m%d")})

# ------------------------------------------------------------------------


class ExporterOldAAA(Exporter):

    """the old AAA format.
    """

    def __init__(self,
                 selector_predictions,
                 selector_genes, selector_mrnas, selector_cds,
                 options,
                 selector_contigs=None,
                 selector_disruptions=None):

        Exporter.__init__(self,
                          selector_predictions,
                          selector_genes, selector_mrnas, selector_cds,
                          options, do_map=False,
                          selector_contigs=selector_contigs,
                          selector_disruptions=selector_disruptions)

        self.mMethod = self.mOptions.method_id

    # ------------------------------------------------------------------------
    def getGeneIdentifier(self, gene_id):
        return "%s_%s_%s" % (self.mOptions.group_id, self.mOptions.species_id, gene_id)

    # ------------------------------------------------------------------------
    def getMRNAIdentifier(self, prediction_id, gene_id):
        if gene_id == "0":
            return "%s.%s" % (prediction_id, prediction_id)
        else:
            return "%s.%s" % (gene_id, prediction_id)

    # ------------------------------------------------------------------------
    def writeHeader(self, outfile):

        outfile.write("""##gff-version   3
#species: %(species)s
#assembly-id: %(assembly)s
#annotation-group-id: %(group_id)s_%(method_id)s
#algorithm: %(group_id)s_%(method_id)s = Oxford gene prediction pipeline based on the program exonerate.
#templates: Gene models from FlyBase 4.2.1 obtained from ENSEMBL version 37
#authors: Andreas Heger and Chris Ponting at firstname.lastname@anat.ox.ac.uk
#date: %(date)s
""" % { 'species': self.mOptions.species,
            'assembly': self.mOptions.assembly_id,
            'group_id': self.mOptions.group_id,
            'method_id': self.mOptions.method_id,
            'date': time.strftime("%Y%m%d")})

# ------------------------------------------------------------------------


class ExporterGBrowser(Exporter):

    def __init__(self,
                 selector_predctions,
                 selector_genes, selector_mrnas, selector_cds,
                 options,
                 with_info=False,
                 selector_contigs=None,
                 selector_disruptions=None):

        Exporter.__init__(self,
                          selector_predictions,
                          selector_genes, selector_mrnas, selector_cds,
                          options, do_map=False,
                          selector_contigs=selector_contigs,
                          selector_disruptions=selector_disruptions)

        self.mWithInfo = True

    # ------------------------------------------------------------------------
    def getGeneIdentifier(self, gene_id):
        return "%s_%s" % (self.mOptions.species_id, gene_id)

    # ------------------------------------------------------------------------
    def getMRNAIdentifier(self, prediction_id, gene_id):
        if gene_id == "0":
            return "%s.%s" % (prediction_id, prediction_id)
        else:
            return "%s.%s" % (gene_id, prediction_id)

    # ------------------------------------------------------------------------
    def getInfoGene(self, gene_id, row):
        """get gene info."""
        return " ; ".join(("Gene %s" % gene_id,))

    # ------------------------------------------------------------------------
    def getInfoMRNA(self, prediction_id, gene_id, row):
        """get mrna info."""
        if self.mWithInfo:
            return " ; ".join(("mRNA %s" % prediction_id,
                               "Gene %s" % gene_id,
                               ":".join(map(str, (row['query_token'],
                                                  row['query_coverage'],
                                                  row['percent_identity'],
                                                  row['nframeshifts'],
                                                  row['nstopcodons']))),
                               "Status %s" % row['class']))
        else:
            return " ; ".join(("mRNA %s" % prediction_id,
                               "Gene %s" % gene_id))

    # ------------------------------------------------------------------------
    def getInfoCDS(self, prediction_id, gene_id, row):
        """get cds info."""
        if self.mWithInfo:
            return " ; ".join(("mRNA %s" % prediction_id,
                               "isOk %s" % row['is_ok'],
                               "pide %i" % row['percent_identity']))
        else:
            return " ; ".join(("mRNA %s" % prediction_id,))

    # ------------------------------------------------------------------------
    def writeHeader(self, outfile):

        outfile.write("""##gff-version   3
#species: %(species)s
#assembly-id: %(assembly)s
#annotation-group-id: %(group_id)s_%(method_id)s
#algorithm: %(group_id)s_%(method_id)s = Oxford gene prediction pipeline based on the program exonerate.
#templates: Gene models from FlyBase 4.2.1 obtained from ENSEMBL version 37
#authors: Andreas Heger and Chris Ponting at firstname.lastname@anat.ox.ac.uk
#date: %(date)s
""" % { 'species': self.mOptions.species,
            'assembly': self.mOptions.assembly_id,
            'group_id': self.mOptions.group_id,
            'method_id': self.mOptions.method_id,
            'date': time.strftime("%Y%m%d")})

# ------------------------------------------------------------------------


class ExporterGTF(Exporter):

    def __init__(self,
                 selector_predctions,
                 selector_genes, selector_mrnas, selector_cds,
                 options,
                 with_info=False,
                 selector_contigs=None,
                 selector_disruptions=None):

        Exporter.__init__(self,
                          selector_predictions,
                          selector_genes, selector_mrnas, selector_cds,
                          options, do_map=False,
                          selector_contigs=selector_contigs,
                          selector_disruptions=selector_disruptions)

        self.mWithInfo = True

        self.mOptions.tracks = ("genes", "mrnas", "cds")

        self.mWriteGeneTrack = False
        self.mWriteMRNATrack = False
        self.mWriteCDSTrack = True
        self.mWriteExonTrack = True
        self.mWriteStopCodon = True
        self.mWriteStartCodon = True

    # ------------------------------------------------------------------------
    def getGeneIdentifier(self, gene_id):
        return "%s_%s" % (self.mOptions.species_id, gene_id)

    # ------------------------------------------------------------------------
    def getMRNAIdentifier(self, prediction_id, gene_id):
        return "%s.%s" % (gene_id, prediction_id)

    # ------------------------------------------------------------------------
    def getInfoCDS(self, prediction_id, gene_id, row):
        """get cds info."""
        return 'gene_id "%s"; transcript_id "%s"; exon_number "%i"' % (gene_id, prediction_id, row['exon_id'])

    # ------------------------------------------------------------------------
    def writeHeader(self, outfile):

        outfile.write("""##gff-version   3
#species: %(species)s
#assembly-id: %(assembly)s
#annotation-group-id: %(group_id)s_%(method_id)s
#algorithm: %(group_id)s_%(method_id)s = Oxford gene prediction pipeline based on the program exonerate.
#templates: Gene models from FlyBase 4.2.1 obtained from ENSEMBL version 37
#authors: Andreas Heger and Chris Ponting at firstname.lastname@anat.ox.ac.uk
#date: %(date)s
""" % { 'species': self.mOptions.species,
            'assembly': self.mOptions.assembly_id,
            'group_id': self.mOptions.group_id,
            'method_id': self.mOptions.method_id,
            'date': time.strftime("%Y%m%d")})


# ------------------------------------------------------------------------
def writeFasta(filenames, options):
    """write fasta files.

    This procedure makes sure that contig names in the index
    correspond to the fasta headers. This is not ensured, because
    the index building allows to extract the contig name from
    anywhere on the fasta header.
    """

    for filename in filenames:
        if filename[-6:] == ".fasta":
            filename = filename[:-6]
        if filename[-3:] == ".fa":
            filename = filename[:-3]
        fasta = IndexedFasta.IndexedFasta(filename)
        contigs = fasta.getContigs()
        for contig in contigs:
            sequence = fasta.getSequence(contig, "+", 0, 0)
            options.stdout.write(">%s\n%s\n" % (contig, sequence))

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/export_predictions.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-s", "--set", dest="set",
                      help="set to use.", type="choice",
                      choices=("raw", "filtered", "clean", "optic", "orthologs", "max-ks", "max-ks-quality", "filter-predictions", "filter-genes"))

    parser.add_option("-f", "--format", dest="format",
                      help="format to use.", type="choice",
                      choices=("oldAAA", "AAA", "GBrowser", "GTF", "fasta"))

    parser.add_option("-b", "--filename-batch", dest="filename_batch", type="string",
                      help="filename with batch information.")
    parser.add_option("--filter-tsv-file", dest="filename_filter", type="string",
                      help="filename with predictions/genes to export.")
    parser.add_option("-o", "--filename-pattern-output", dest="filename_pattern_output", type="string",
                      help="filename with one %s for batch output.")
    parser.add_option("-m", "--filename-pattern-map", dest="filename_pattern_map", type="string",
                      help="filename with one %s for maps.")
    parser.add_option("-i", "--first-id", dest="first_id", type="int",
                      help="first number to be used for identifiers.")
    parser.add_option("--tablename-orthologs", dest="tablename_orthologs", type="string",
                      help="tablename with orthology information. Only orthologs are output.")
    parser.add_option("--schema", dest="schema", type="string",
                      help="schema to use.")
    parser.add_option("--species", dest="species", type="string",
                      help="full species name.")
    parser.add_option("--species-id", dest="species_id", type="string",
                      help="species abbreviation.")
    parser.add_option("--group-id", dest="group_id", type="string",
                      help="group abbreviation.")
    parser.add_option("--method-id", dest="method_id", type="string",
                      help="method id to use.")
    parser.add_option("--tracks", dest="tracks", type="string",
                      help="tracks to export [contigs|genes|mrnas|cds].")
    parser.add_option("--max-ks", dest="max_ks", type="float",
                      help="maximum ks [0]. Use for set max-ks.")

    parser.set_defaults(
        set="clean",
        format="AAA",
        filename_batch=None,
        filename_pattern_output="%s.gff",
        filename_pattern_map=None,
        filename_filter=None,
        species="Drosophila pseudoobscura",
        assembly_id="dpse_caf1",
        abbreviation="dpse",
        species_id="GA",
        method_id="GPI",
        group_id="OXFD",
        schema="dpse_vs_dmel8",
        tablename_predictions="predictions",
        tablename_genes="genes",
        tablename_redundant="redundant",
        tablename_exons="exons",
        tablename_contigs="contigs",
        tablename_quality="quality",
        tablename_orthologs=None,
        tablename_kaks="kaks",
        tablename_disruptions="disruptions",
        first_id=2000000,
        max_ks=0,
        quality_clean_set='CG,PG,SG',
        quality_optic_set='CG,SG,PG,RG,CP,SP,PP',
        tracks="contigs,genes,mrnas,cds",
    )

    (options, args) = E.Start(parser,
                              add_database_options=True,
                              add_pipe_options=True)

    options.quality_clean_set = options.quality_clean_set.split(",")
    options.quality_optic_set = options.quality_optic_set.split(",")
    options.tracks = options.tracks.split(",")

    dbhandle = pgdb.connect(options.psql_connection)

    # pre-processing: set filtering options
    selector_genes = SelectorGenes(dbhandle, options)
    selector_mrnas = SelectorMRNA(dbhandle, options)
    selector_cds = SelectorCDS(dbhandle, options)

    if options.set == "raw":
        selector_predictions = SelectorPredictionsRaw(dbhandle, options)
    elif options.set == "filtered":
        selector_predictions = SelectorPredictionsFiltered(dbhandle, options)
    elif options.set == "clean":
        selector_predictions = SelectorPredictionsClean(dbhandle, options)
    elif options.set == "optic":
        selector_predictions = SelectorPredictionsOptic(dbhandle, options)
    elif options.set == "orthologs":
        selector_predictions = SelectorPredictionsOrthologs(dbhandle, options)
    elif options.set == "max-ks":
        selector_predictions = SelectorPredictionsMaxKs(dbhandle, options)
    elif options.set == "max-ks-quality":
        selector_predictions = SelectorPredictionsMaxKsQuality(
            dbhandle, options)
    elif options.set == "filter-predictions":
        selector_predictions = SelectorPredictionsFilterPredictions(
            dbhandle, options)
    elif options.set == "filter-genes":
        selector_predictions = SelectorPredictionsFilterGenes(
            dbhandle, options)

    selector_disruptions = SelectorDisruptions(dbhandle, options)

    if options.format == "AAA":
        exporter_class = ExporterAAA
        selector_contigs = None
    elif options.format == "oldAAA":
        exporter_class = ExporterOldAAA
        selector_contigs = None
    elif options.format == "GBrowser":
        exporter_class = ExporterGBrowser
        selector_contigs = SelectorContigs(dbhandle, options)
    elif options.format == "GTF":
        exporter_class = ExporterGTF
        selector_contigs = None

    elif options.format == "fasta":

        writeFasta(args, options)
        E.Stop()
        sys.exit(0)

    if options.filename_batch:
        infile = open(options.filename_batch)
        for line in infile:
            if line[0] == "#":
                continue
            options.species, options.abbreviation, options.species_id, options.assembly_id, options.schema = line[
                :-1].split("\t")[:5]
            filename = options.filename_pattern_output % options.assembly_id
            outfile = open(filename, "w")

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# writing data from %s to %s\n" % (options.schema, filename))
                options.stdlog.flush()

            exporter = exporter_class(selector_predictions,
                                      selector_genes, selector_mrnas, selector_cds, options,
                                      selector_contigs=selector_contigs,
                                      selector_disruptions=selector_disruptions)
            exporter.setMethod(options.method_id)
            map_gene_id, map_prediction_id = exporter.writeGFF(
                outfile, dbhandle, options)

            outfile.close()

            if options.filename_pattern_map:
                exporter.writeMap(options, "genes", map_gene_id)
                exporter.writeMap(options, "predictions", map_prediction_id)

        infile.close()

    else:
        outfile = options.stdout

        exporter = exporter_class(selector_predictions,
                                  selector_genes, selector_mrnas, selector_cds, options,
                                  selector_contigs=selector_contigs,
                                  selector_disruptions=selector_disruptions)

        exporter.setMethod(options.method_id)

        map_gene_id, map_prediction_id = exporter.writeGFF(outfile)

        if options.filename_pattern_map and map_gene_id and map_prediction_id:
            WriteMap(options, "genes", map_gene_id)
            WriteMap(options, "predictions", map_prediction_id)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
