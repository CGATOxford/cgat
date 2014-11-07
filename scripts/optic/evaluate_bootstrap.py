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
optic/evaluate_bootstrap.py - evaluate bootstrap results against reference
====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

read in a tree in newick format and compare bootstrap as given
by results if they agree with this tree.

The input format is the output from phylip's >consense< program.

Options::

    -h, --help                      print this message.
    -v, --verbose=                  loglevel.
    -p, --pattern-species=          regex pattern to extract species from identifier
    -f, --file-tree=                filename with reference-tree
    -r, --reference-tree            reference tree in the command line (no spaces!)
    -c, --min-cluster-support       cluster according to species tree (minimum support)
    -e, --min-report-support        report consistent/inconsistent clades (minimum support). Set
                                    to value below 0, if you want to process lists until first
                                    inconsistent edge.
    -c, --file-clusters=            filename with clusters
    -i, --file-inconsistencies=     filename with inconsistencies
    -s, --file-subtrees=            filename with subtrees
    -x, --column-prefix                    cluster prefix to use

Usage
-----

Example::

   python optic/evaluate_bootstrap.py --help

Type::

   python optic/evaluate_bootstrap.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import time
import sets

from Bio.Nexus import Nexus
from Bio.Nexus.Nodes import Node

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

param_loglevel = 1

param_long_options = ["verbose=", "help", "pattern-species=",
                      "reference-tree=", "file-tree=", "cluster",
                      "min-cluster-support=", "min-report-support=",
                      "file-clusters=", "file-inconsistencies=",
                      "file-subtrees=",
                      "prefix=", "version"]

param_short_options = "v:hp:f:r:c:i:x:s:"

param_pattern_species = "^([^@|:]+)[@|:]"

param_reference_tree = None

param_cluster = False

param_min_cluster_support = 90
param_min_report_support = 80

param_filename_clusters = None
param_filename_inconsistencies = None
param_filename_subtrees = None

param_prefix_cluster = "0"


class Results:
    ntests = 0
    norgs = 0
    notus = 0
    ntrivial = 0
    nfull = 0
    nskipped = 0
    nconsistent_external = 0
    ninconsistent_external = 0
    nconsistent_internal = 0
    ninconsistent_internal = 0
    nfixed_consistent_external = 0
    nfixed_inconsistent_external = 0
    nfixed_consistent_internal = 0
    nfixed_inconsistent_internal = 0
    ntotal_external = 0
    ntotal_internal = 0
    mask = None
    mask_id = 0
    inconsistencies = None
    is_fixed = False

    def __init__(self,
                 mask=None,
                 notus=0,
                 norgs=0,
                 mask_id=0):

        self.notus = notus

        self.setOrgs(norgs)

        if mask:
            self.mask = mask
        else:
            self.mask = [0] * self.notus

        self.mask_id = mask_id
        self.inconsistencies = []
        self.subtrees = {}

    def printHeader(self):
        """print header to go with the str() function."""
        return string.join(("mask", "mask_id",
                            "norgs", "notus", "ntests",
                            "nci", "nii", "nti",
                            "nce", "nie", "nte",
                            "nfci", "nfii",
                            "nfce", "nfie",
                            "nfull", "ntriv", "nskip"), "\t")

    def __str__(self):
        return string.join(map(str, (
            string.join(map(str, self.mask), ""),
            self.mask_id,
            self.norgs, self.notus, self.ntests,
            self.nconsistent_internal, self.ninconsistent_internal, self.ntotal_internal,
            self.nconsistent_external, self.ninconsistent_external, self.ntotal_external,
            self.nfixed_consistent_internal, self.nfixed_inconsistent_internal,
            self.nfixed_consistent_external, self.nfixed_inconsistent_external,
            self.nfull, self.ntrivial, self.nskipped)), "\t")

    def setOrgs(self, norgs):

        self.norgs = norgs
        self.ntotal_external = self.norgs
        self.ntotal_internal = max(0, 2 * (self.norgs - 3))

    def isOk(self):
        """returns True, if a good subset."""
        return self.ninconsistent_internal == 0 and self.ninconsistent_external == 0

    def isOverlap(self, other):
        return reduce(lambda x, y: x + y, map(lambda x, y: x * y, self.mask, other.mask)) > 0

    def getSize(self):
        return reduce(lambda x, y: x + y, self.mask)

    def getLength(self):
        return len(self.mask)

    def printSummary(self, format="tab"):
        """return summary string."""

        t1 = max(self.ntests, 1)
        t2 = max(self.ntotal_internal + self.ntotal_external, 1)

        if format == "pretty":
            return string.join(map(str, (
                "norgs=%i" % self.norgs,
                "notus=%i" % self.notus,
                "ptest=%5.4f" %
                (float(
                    self.nconsistent_internal + self.nconsistent_external) / t1),
                "ptotal=%5.4f" %
                (float(
                    self.nconsistent_internal + self.nconsistent_external) / t2),
                "ftotal=%5.4f" %
                (float(self.nfixed_consistent_internal +
                 self.nfixed_consistent_external) / t2),
            )), "\t")
        elif format == "tab":
            return string.join(map(str, (
                "%i" % self.norgs,
                "%i" % self.notus,
                "%5.4f" %
                (float(
                    self.nconsistent_internal + self.nconsistent_external) / t1),
                "%5.4f" %
                (float(
                    self.nconsistent_internal + self.nconsistent_external) / t2),
                "%5.4f" %
                (float(self.nfixed_consistent_internal +
                 self.nfixed_consistent_external) / t2),
            )), "\t")
        else:
            raise "Unknown format %s" % format

    def fix(self):
        """fix current values of consistencies and inconsistencies."""
        if self.is_fixed:
            return
        self.is_fixed = True
        self.nfixed_consistent_internal, self.nfixed_inconsistent_internal = \
            self.nconsistent_internal, self.ninconsistent_internal

        self.nfixed_consistent_external, self.nfixed_inconsistent_external = \
            self.nconsistent_external, self.ninconsistent_external

    def increment_external(self, is_within_support, is_consistent, present):

        if is_within_support:
            self.ntests += 1

            if is_consistent:
                self.nconsistent_external += 1
            else:
                self.ninconsistent_external += 1
                self.inconsistencies.append(present)

        if not self.is_fixed:

            if is_consistent:
                self.nfixed_consistent_external += 1
            else:
                self.nfixed_inconsistent_external += 1
                self.is_fixed = True

    def increment_internal(self, is_within_support, is_consistent, present):

        if is_within_support:
            self.ntests += 1

            if is_consistent:
                self.nconsistent_internal += 1
            else:
                self.ninconsistent_internal += 1
                self.inconsistencies.append(present)

        if not self.is_fixed:

            if is_consistent:
                self.nfixed_consistent_internal += 1
            else:
                self.nfixed_inconsistent_internal += 1
                self.is_fixed = True

    def addSubtree(self, present):
        """add subtree to list of subtrees."""

        p = string.join(map(str, present), "")
        if p not in self.subtrees:
            self.subtrees[p] = 0

        self.subtrees[p] += 1

##########################################################################


def Prune(tree, taxon):
    """Prunes a terminal taxon from the tree.

    id_of_previous_node = prune(tree,taxon)
    If taxon is from a bifurcation, the connecting node will be collapsed
    and its branchlength added to remaining terminal node. This might be no
    longer a meaningful value'
    """

    id = tree.search_taxon(taxon)
    if id is None:
        raise TreeError('Taxon not found: %s' % taxon)
    elif id not in tree.get_terminals():
        raise TreeError('Not a terminal taxon: %s' % taxon)
    else:
        prev = tree.unlink(id)
        tree.kill(id)
        if not prev == tree.root and len(tree.node(prev).succ) == 1:
            succ = tree.node(prev).get_succ()
            tree.collapse(prev)

        return prev

##########################################################################


def GetPrunedReferenceTree(mask, present_orgs, reference_tree):

    # reread and process species tree
    # has to be done for every new pass, because
    # the tree is modified later on (and I haven't found
    # a copy mechanism (because I did not look)).
    nexus = TreeTools.Newick2Nexus(reference_tree)
    reference_tree = nexus.trees[0]

    ###########################################################################
    # prune reference tree and keep only those taxa, which are present in the
    # cluster.
    for nx in reference_tree.get_terminals():
        otu = reference_tree.node(nx).get_data().taxon
        if otu not in present_orgs:
            Prune(reference_tree, otu)

    if param_loglevel >= 3:
        print "# pruned reference tree for %s:" % (",".join(present_orgs.keys()))
        reference_tree.display()

    return reference_tree

##########################################################################


def IsMonophyletic(reference_tree, t1, t2):
    """decide whether the partition given by t1 and t2 are monophyletic.
    """

    is_overlap = False
    for nx in reference_tree.get_terminals():
        t = reference_tree.node(nx).get_data().taxon
        if t in t1 and t in t2:
            is_overlap = True
            break
    if is_overlap:
        is_monophyletic1 = False
        is_monophyletic2 = False
    else:
        reference_tree.root_with_outgroup(t2)
        is_monophyletic1 = reference_tree.is_monophyletic(t1) >= 0

        reference_tree.root_with_outgroup(t1)
        is_monophyletic2 = reference_tree.is_monophyletic(t2) >= 0

    # root with something not in t1
    for nx in reference_tree.get_terminals():
        otu = reference_tree.node(nx).get_data().taxon
        if otu not in t1:
            reference_tree.root_with_outgroup(otu)
            is_subtree1 = reference_tree.is_monophyletic(t1) >= 0
            break
    else:
        is_subtree1 = False

    # root with something not in t2
    for nx in reference_tree.get_terminals():
        otu = reference_tree.node(nx).get_data().taxon
        if otu not in t2:
            reference_tree.root_with_outgroup(otu)
            is_subtree2 = reference_tree.is_monophyletic(t2) >= 0
            break
    else:
        is_subtree2 = False

    return is_monophyletic1, is_monophyletic2, is_subtree1, is_subtree2

##########################################################################


def AddResults(mask,
               support, min_report_support,
               pattern,
               is_monophyletic, is_subtree,
               taxa,
               present):
    """add evaluation results.

    returns false, in case the pattern is inconsistent with the species tree.
    """

    do_print = True
    is_ok = True

    # do consistency-counting
    if len(taxa) == 1:

        mask.ntrivial += 1
        do_print = False

    elif len(taxa) == mask.norgs:

        mask.nfull += 1
        do_print = False

    elif len(taxa) == mask.norgs - 1:
        mask.increment_external(support >= min_report_support,
                                is_monophyletic,
                                present)
        is_ok = is_monophyletic
    else:
        mask.increment_internal(support >= min_report_support,
                                is_monophyletic,
                                present)
        is_ok = is_monophyletic

        if is_subtree:
            mask.addSubtree(present)

    if do_print and param_loglevel >= 2:
        print "# " + string.join(map(str, (
            pattern, support,
            is_monophyletic, is_subtree,
            string.join(map(str, present), ""))), "\t")

    return is_ok


def GetOrgs(map_id2org):
    """return hash of present organisms in map."""

    present_orgs = {}
    for x in range(len(map_id2org)):
        org, name, nid = map_id2org[x]
        present_orgs[org] = 1

    return present_orgs


##########################################################################
def SelectMasks(masks, patterns, norgs, map_id2org, min_report_support=90):
    """from a set of evaluated masks, select a best subset.
    """

    new_masks = []

    # Sort by size of mask and start with the smallest
    masks.sort(lambda x, y: cmp(x.getSize(), y.getSize()))

    # parse masks and select those that have no inconsistencies and are not
    # overlapping.
    for mask in masks:
        if mask.isOk():
            is_overlap = False
            for new_mask in new_masks:
                if mask.isOverlap(new_mask):
                    is_overlap = True
                    break
            if not is_overlap:
                new_masks.append(mask)

    # parse masks and select those that have inconsistencies and are not
    # overlapping.
    for mask in masks:
        if not mask.isOk():
            is_overlap = False
            for new_mask in new_masks:
                if mask.isOverlap(new_mask):
                    is_overlap = True
                    break
            if not is_overlap:
                new_masks.append(mask)

    # add a final mask for non-overlapping patterns that covers all
    # remaining otus.
    cover_mask = Results([1] * masks[0].getLength(), mask_id=0)

    for m in new_masks:
        for x in range(len(m.mask)):
            cover_mask.mask[x] = cover_mask.mask[x] & (not m.mask[x])

    if cover_mask.getSize() > 0:

        reference_tree = GetPrunedReferenceTree(
            cover_mask, GetOrgs(map_id2org), param_reference_tree)
        cover_mask.setOrgs(norgs)

        AnalyseMask(cover_mask, patterns, norgs, reference_tree,
                    map_id2org, min_report_support)
        new_masks.append(cover_mask)

    return new_masks

##########################################################################


def AnalyseMask(mask, patterns, norgs, reference_tree, map_id2org, min_report_support=90):
    """analyse partitions with mask."""

    tested = {}

    for support, pattern in patterns:

        t1, t2 = {}, {}
        p1, p2 = [0] * norgs, [0] * norgs

        if param_loglevel >= 4:
            print "# pattern=", pattern

        for x in range(len(pattern)):

            org, name, nid = map_id2org[x]

            if org == "unknown":
                continue

            if mask.mask[x]:
                if pattern[x] == "*":
                    t1[org] = 1
                    p1[nid] = 1
                    if param_loglevel >= 4:
                        print "# *", name, nid
                else:
                    t2[org] = 1
                    p2[nid] = 1
                    if param_loglevel >= 4:
                        print "# .", name, nid

        t1 = t1.keys()
        t2 = t2.keys()
        t1.sort()
        t2.sort()

        # do not test the same partitions twice
        key1 = string.join(map(str, p1), "") + string.join(map(str, p2), "")
        key2 = string.join(map(str, p2), "") + string.join(map(str, p1), "")
        if key1 in tested or key2 in tested:
            continue

        tested[key1] = 1
        tested[key2] = 1

        if len(t1) == 0 or len(t2) == 0:
            continue
        is_monophyletic1, is_subtree1, is_monophyletic2, is_subtree2 = IsMonophyletic(
            reference_tree, t1, t2)

        AddResults(mask, support, min_report_support, pattern,
                   is_monophyletic1, is_subtree1, t1, p1)
        AddResults(mask, support, min_report_support, pattern,
                   is_monophyletic2, is_subtree2, t2, p2)

##########################################################################


def AnalysePatterns(patterns,
                    map_id2org,
                    min_cluster_support=100,
                    min_report_support=90):
    """analyse partitions by comparing to reference tree.

    Prints out for each partition, whether left/right is consistent
    with reference tree or not.

    If there are full complements on either side, print suggested split.

    Prints summary statistics:
    for each consistent partition:
        print counts
    """

    # reread and process species tree
    # has to be done for every new pass, because
    # the tree is modified later on (and I haven't found
    # a copy mechanism (because I did not look)).
    nexus = TreeTools.Newick2Nexus(param_reference_tree)
    reference_tree = nexus.trees[0]

    norgs = len(reference_tree.get_terminals())
    notus = len(patterns[0][1])

    # complement patterns with single species patterns:
    patterns.reverse()
    for x in range(notus):
        pattern = ["."] * notus
        pattern[x] = "*"
        patterns.append((100, string.join(pattern, "")))
    patterns.reverse()

    ##########################################################################
    # first pass: separate well supported full species trees
    masks = []
    present_orgs = {}
    mask_id = 0

    for support, pattern in patterns:

        t1, t2, i1, i2 = {}, {}, [], []

        for x in range(len(pattern)):
            org, name, nid = map_id2org[x]
            if org == "unknown":
                continue
            present_orgs[org] = 1
            if pattern[x] == "*":
                t1[org] = 1
                i1.append(name)
            else:
                t2[org] = 1
                i2.append(name)

        t1 = t1.keys()
        t2 = t2.keys()
        t1.sort()
        t2.sort()

        if param_loglevel >= 4:
            print "# ", pattern, len(t1), len(t2), i1, i2
            sys.stdout.flush()

        if len(t1) == len(t2) and \
                len(t1) == norgs and \
                support >= min_cluster_support:

            mask1, notus1 = [], 0
            mask2, notus2 = [], 0

            for x in range(len(pattern)):
                if pattern[x] == "*":
                    notus1 += 1
                    mask1.append(1)
                    mask2.append(0)
                else:
                    notus2 += 1
                    mask1.append(0)
                    mask2.append(1)

            mask_id += 1
            masks.append(Results(mask1, notus1, len(t1), mask_id=mask_id))
            mask_id += 1
            masks.append(Results(mask2, notus2, len(t2), mask_id=mask_id))

            if param_loglevel >= 2:
                print "# split\tfull\t%i\t%s\t%i\t%i\t%s" % (support, string.join(map(str, mask1), ""), notus1, len(t1), string.join(i1, ";"))
                print "# split\tfull\t%i\t%s\t%i\t%i\t%s" % (support, string.join(map(str, mask2), ""), notus2, len(t2), string.join(i2, ";"))

    # add full mask
    if len(masks) == 0:
        masks.append(Results([1] * notus, notus, len(present_orgs), mask_id=1))

    ##########################################################################
    # second pass: check subtrees for each mask
    # external: edges leading to external nodes (i.e., leaves): total number = norgs
    # internal: all other edges: maximum number = 2 * (2 * norgs - 3 - norgs) = 2 * (norgs - 3)
    # 1st factor 2: two directions
    # 2nd factor: 2n-3 is number of edges in unrooted tree.
    # 3rd factor: -n = number of external edges
    for mask in masks:
        reference_tree = GetPrunedReferenceTree(
            mask, GetOrgs(map_id2org), param_reference_tree)
        AnalyseMask(
            mask, patterns, norgs, reference_tree, map_id2org, min_report_support)

    if param_loglevel >= 1:
        print "# partitions after evaluation:"
        print "#", Results().printHeader()
        for m in masks:
            print "#", str(m)

    reference_tree = GetPrunedReferenceTree(
        mask, GetOrgs(map_id2org), param_reference_tree)
    new_masks = SelectMasks(masks, patterns,
                            norgs, map_id2org,
                            min_report_support)

    if param_loglevel >= 1:
        print "# partitions after selection:"
        print "#", Results().printHeader()
        for m in new_masks:
            print "#", str(m)

    return new_masks

##########################################################################


def WriteClusters(outfile, mask, map_id2org, cluster_id):
    """write new cluster.
    """

    outfile.write(">cluster# %s_%s\n" % (cluster_id, mask.mask_id))
    for x in range(len(mask.mask)):
        if mask.mask[x]:
            outfile.write("%s\t%s_%s\n" %
                          (map_id2org[x][1], cluster_id, mask.mask_id))

##########################################################################


def WriteInconsistencies(outfile, mask, map_id2org, cluster_id):
    """write inconsistencies per cluster. Do not write empty entries and
    do not write the covering mask.
    """

    if len(mask.inconsistencies) == 0 or mask.mask_id == 0:
        return

    outfile.write(">cluster# %s_%s\n" % (cluster_id, mask.mask_id))
    mask.inconsistencies.sort()

    for p in mask.inconsistencies:
        outfile.write("%s\t%s_%s\n" %
                      (string.join(map(str, p), ""), cluster_id, mask.mask_id))


##########################################################################
def WriteSubtrees(outfile, mask, map_id2org, cluster_id):
    """write inconsistencies per cluster. Do not write empty entries and
    do not write the covering mask.
    """

    if len(mask.subtrees) == 0:
        return

    outfile.write(">cluster# %s_%s\n" % (cluster_id, mask.mask_id))
    keys = mask.subtrees.keys()
    keys.sort()

    for k in keys:
        outfile.write("%s\t%s_%s\t%i\n" %
                      (k, cluster_id, mask.mask_id, mask.subtrees[k]))

##########################################################################


def ParseTree(reference_tree, rx_species):

    nexus = TreeTools.Newick2Nexus(reference_tree)
    reference_tree = nexus.trees[0]
    if param_loglevel >= 3:
        print "# reference tree:"
        reference_tree.display()

    map_taxon2id = {}
    for nx in reference_tree.get_terminals():
        otu = reference_tree.node(nx).get_data().taxon
        map_taxon2id[otu] = len(map_taxon2id)
        if param_loglevel >= 2:
            print "# %s\t%i" % (otu, map_taxon2id[otu])
    map_taxon2id["unknown"] = len(map_taxon2id)

    return reference_tree, map_taxon2id

##########################################################################


def ReadMapId2Org(infile, rx_species, map_taxon2id):
    """read map section from stdin.
    """

    map_id2org = []

    infile.readline()

    x = 0
    map_otu2id = {}
    while 1:

        line = infile.readline()
        if not line:
            break

        try:
            id, otu = re.search("(\d)+.\s+(\S+)", line[:-1]).groups()
        except AttributeError:
            break

        try:
            org = rx_species.search(otu).groups()[0]
        except AttributeError:
            org = "unknown"

        map_otu2id[otu] = len(map_id2org)
        map_id2org.append((org, otu, map_taxon2id[org]))

        x += 1
        if x == 10:
            x = 0
            infile.readline()

    if param_loglevel >= 2:
        print "# map_id2org=", map_id2org
        sys.stdout.flush()

    return map_id2org, map_otu2id


##########################################################################
def ReadPatterns(infile):
    """read patterns section from stdin.
    """
    line = infile.readline()

    line = infile.readline()
    num_samples = int(re.search("(\d+)", line[:-1]).groups()[0])

    line = infile.readline()

    if param_loglevel >= 2:
        print "# num_samples=", num_samples
        sys.stdout.flush()

    patterns = []

    while 1:
        line = infile.readline()
        if not line:
            break

        data = re.split("\s+", line[:-1])
        if len(data) < 2:
            break

        try:
            pattern = string.join(data[:-1], "")
            support = int(float(data[-1]))
        except ValueError:
            continue

        if param_loglevel >= 2:
            print "# pattern\t%s\tsupport\t%i" % (pattern, support)

        patterns.append((support, pattern))

    return patterns, num_samples

##########################################################################
##########################################################################
##########################################################################


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ("-p", "--pattern-species"):
            param_pattern_species = a
        elif o in ("-t", "--reference-tree"):
            param_reference_tree = a
        elif o in ("-f", "--file-tree"):
            param_reference_tree = filter(
                lambda x: x[0] != "#", open(a, "r").readlines())
        elif o in ("-m", "--min-cluster-support"):
            param_min_cluster_support = int(a)
        elif o in ("-e", "--min-report-support"):
            param_min_report_support = int(a)
        elif o in ("-c", "--file-clusters"):
            param_filename_clusters = a
        elif o in ("-i", "--file-inconsistencies"):
            param_filename_inconsistencies = a
        elif o in ("-s", "--file-subtrees"):
            param_filename_subtrees = a
        elif o in ("-x", "--column-prefix"):
            param_prefix_cluster = a

    print E.GetHeader()
    print E.GetParams()

    rx_species = re.compile(param_pattern_species)

    # read and process reference tree
    if not param_reference_tree:
        print USAGE
        print "please supply reference tree"
        sys.exit(1)

    reference_tree, map_taxon2id = ParseTree(param_reference_tree, rx_species)

    # open output files
    if param_filename_clusters:
        outfile_clusters = open(param_filename_clusters, "a")
    else:
        outfile_clusters = sys.stdout

    if param_filename_inconsistencies:
        outfile_inconsistencies = open(param_filename_inconsistencies, "a")
    else:
        outfile_inconsistencies = sys.stdout

    if param_filename_subtrees:
        outfile_subtrees = open(param_filename_subtrees, "a")
    else:
        outfile_subtrees = sys.stdout

    cluster_id = param_prefix_cluster

    print """# Fields for tabbed cluster output:
# part_id:              partition id  
# notus:                number of otus in cluster
# norgs:                number of orgs in cluster
# ptest:                ratio of congruent partitions / tested partitions
# ptotal:               ratio of congruent partitions / all partitions
# ftotal:               ratio of fixed congrued partitions / all partitions"""

    print """# Fields for cluster summary output:
# cluster_id:           cluster id
# nmaps:                number of mapped ids
# nsamples:             number of samples
# nsupp:                number of partitions above support
# nparts:               number of partitions
# nmasks:               number of masks
"""

    while 1:
        line = sys.stdin.readline()
        if not line:
            break
        if line[0] == "#":
            continue

        # new cluster entry
        if line[0] == ">":

            cluster_id = re.search(">cluster# (\S+)", line[:-1]).groups()[0]

        elif re.match("Species in order:", line):

            map_id2org, map_otu2id = ReadMapId2Org(
                sys.stdin, rx_species, map_taxon2id)
            continue

        elif re.match("Sets included in the consensus tree", line):

            patterns, num_samples = ReadPatterns(sys.stdin)

            if len(patterns) == 0:
                print "# no patterns, analysis skipped"
                continue

            masks = AnalysePatterns(patterns,
                                    map_id2org,
                                    param_min_cluster_support,
                                    param_min_report_support)

            for mask in masks:
                WriteClusters(outfile_clusters, mask, map_id2org, cluster_id)

            for mask in masks:
                WriteInconsistencies(
                    outfile_inconsistencies, mask, map_id2org, cluster_id)

            for mask in masks:
                WriteSubtrees(outfile_subtrees, mask, map_id2org, cluster_id)

            for mask in masks:
                print "%s_%s" % (cluster_id, mask.mask_id) + "\t" + mask.printSummary()

            print string.join(map(str, ("summary", cluster_id, len(map_id2org), num_samples,
                                        len(filter(
                                            lambda x: x[0] >= param_min_report_support, patterns)),
                                        len(patterns), len(masks))), "\t")

    if param_filename_clusters:
        outfile_clusters.close()

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
