'''TreeReconcilation.py -

These methods were taken from analyze_orthology_multiple.py
as shared with analyze_genetrees.py
'''

from types import *

import CGAT.TreeTools as TreeTools
import numpy


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def readLocations(infile, extract_species):
    """read map of ids to locations from file.

    Assign mappings both in terms of transcripts and genes.
    """

    map_id2location = {}
    for line in infile:
        if line[0] == "#":
            continue
        l = Location()
        data = line[:-1].split("\t")[:5]
        l.contig, l.strand, l.mFrom, l.mTo = \
            data[1], data[2], int(data[3]), int(data[4])

        # translate chromosomal names like dsim_chrU to chrU
        # such that filtering for junk chromosomes works.
        # The pattern here works for the flies.
        # NB: I had not translated all names yet, so this is a patch
        # example: dyak: chr2L, but dsim: dsim_chr2L
        fragments = data[1].split("_")
        if fragments[0][0] == "d" and len(fragments[0]) == 4:
            # assume first part is part of species name and will be removed
            l.mShortName = "_".join(fragments[1:])
        else:
            l.mShortName = data[1]

        map_id2location[data[0]] = l

    # sort locations and number them to easily find adjacent genes
    keys = map_id2location.keys()

    # first sort keys by species (this assumes, of course,
    # that identifiers start with a species name)
    keys.sort()

    def processKeys(keys):
        # just one line to sort, but very inefficient
        keys.sort(lambda x, y: cmp((map_id2location[x].contig,
                                    map_id2location[x].mFrom),
                                   (map_id2location[y].contig,
                                    map_id2location[y].mFrom)))
        last_contig = None
        last_to = 0
        x = 0
        for k in keys:
            l = map_id2location[k]
            if last_contig != l.contig or last_to < l.mFrom:
                x += 1
                last_to = l.mTo
            else:
                last_to = max(last_to, l.mTo)

            last_contig = l.contig
            l.mId = x

    first = 0
    last_species = None
    for k in range(len(keys)):
        this_species = extract_species(keys[k])
        if last_species != this_species:
            processKeys(keys[first:k])
            first = k
        last_species = this_species

    processKeys(keys[first:len(keys)])

    return map_id2location


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def filterTree(tree, options, map_id2location=None):
    """apply location and type filter to tree.

    if outgroups are defined, they are not removed.
    """

    otus = TreeTools.GetTaxa(tree)

    to_remove = set()
    if options.remove_unplaced:
        tt = set()
        for id in otus:
            if id not in map_id2location:
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# WARNING: unknown location for id %s.\n" % id)
                continue

            if map_id2location[id].mShortName.lower() in MAP_CONTIG2JUNK:
                to_remove.add(id)
                tt.add(id)

        if options.loglevel >= 3:
            options.stdlog.write("# tree %s: removing %i entries because of location: %s\n" %
                                 (tree.name, len(tt), ";".join(tt)))

    new_otus = list(set(otus).difference(to_remove))

    if len(new_otus) != len(otus):

        TreeTools.PruneTree(tree, new_otus, keep_distance_to_root=True)

    if options.loglevel >= 1:
        options.stdlog.write("# tree %s: filtering: before=%i, remove=%i, after=%i, final=%i\n" %
                             (tree.name, len(otus), len(to_remove), len(new_otus), len(TreeTools.GetTaxa(tree))))
        options.stdlog.flush()

##         quality_codes = MAP_FILTER2QUALITY[options.filter_quality]
##         remove_otus1 = filter( lambda x: x.split(options.separator)[3] not in quality_codes, otus )

# if options.loglevel >= 3:
# options.stdlog.write("# tree %s: removing because of quality: %s\n" % \
# (tree.name, ";".join(remove_otus1)))

##         remove_otus = remove_otus.union(remove_otus1)

# remove_otus2=set()
# if options.filter_location != "all":
##         remove_otus2 = []
# for id in otus:
# if id not in map_id2location:
# if options.loglevel >= 1:
# options.stdlog.write("# WARNING: unknown location for id %s.\n" % id )
# continue
# if map_id2location[id] not in MAP_CONTIG2JUNK:
##                 remove_otus2.append( id )

# if options.loglevel >= 3:
# options.stdlog.write("# tree %s: removing because of location: %s\n" % \
# (tree.name, ";".join(remove_otus2)))

##         remove_otus = remove_otus.union(remove_otus2)

# removing outgroup species from removal list
# if options.outgroup_species:

# find monophyletic trees of outgroup_species
##         outgroup_taxa = filter( lambda x: extract_species(x) in options.outgroup_species, otus)
##         remove_otus = remove_otus.difference( outgroup_taxa )

##     new_otus = list(set(otus).difference(remove_otus))

# if len(new_otus) != len(otus):
##         TreeTools.PruneTree( tree, new_otus )

# if options.loglevel >= 1:
# options.stdlog.write("# tree %s: filtering: before=%i, quality=%i, location=%i, new=%i, final=%i\n" % \
# (tree.name, len(otus), len(remove_otus1), len(remove_otus2), len(new_otus), len(TreeTools.GetTaxa(tree))) )


def rerootTree(gene_tree, extract_species, options):

    otus = TreeTools.GetTaxa(gene_tree)

    # find monophyletic trees of outgroup_species
    try:
        outgroup_taxa = filter(
            lambda x: extract_species(x) in options.outgroup_species, otus)
    except AttributeError:
        raise "error while rerooting tree in tree %s with %s" % (
            gene_tree.name, str(otus))

    if gene_tree.is_monophyletic(outgroup_taxa):
        r = outgroup_taxa
    else:
        r = [outgroup_taxa[0], ]

    if r:
        if options.loglevel >= 1:
            options.stdlog.write("# tree %s: rerooting with %i outgroups:  %s.\n" % (
                gene_tree.name, len(r), ",".join(r)))
            options.stdlog.flush()
    else:
        if options.loglevel >= 1:
            options.stdlog.write(
                "# tree %s: no outgroup found, tree will not be rerooted.\n" % gene_tree.name)
            options.stdlog.flush()

    gene_tree.root_with_outgroup(r)

    if options.loglevel >= 5:
        gene_tree.display()


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def extractSubtrees(tree, extract_species, options):
    """extract subtrees from tree.

    Splits a rooted tree at outgroups or at out-paralogs.

    Returns a list of clusters with its members belonging to each
    subtree.
    """

    nin = [0] * TreeTools.GetSize(tree)
    nout = [0] * TreeTools.GetSize(tree)
    nstop = [False] * TreeTools.GetSize(tree)
    taxa = [set() for x in range(TreeTools.GetSize(tree))]
    otus = [set() for x in range(TreeTools.GetSize(tree))]

    clusters = []

    def update_groups(node_id):
        node = tree.node(node_id)

        if node.succ == []:
            taxa[node_id] = set((extract_species(node.data.taxon),))
            otus[node_id] = set((node.data.taxon,))
            if extract_species(node.data.taxon) in options.outgroup_species:
                nout[node_id] = 1
            else:
                nin[node_id] = 1
        else:
            a, b = node.succ
            oa = nout[a] > 0
            ob = nout[b] > 0
            ia = nin[a] > 0
            ib = nin[b] > 0
            overlap = len(taxa[a].intersection(taxa[b])) > 0

            # merge, if
            # * both have outgroups but no ingroups
            # * either have ingroups but no outgroups.
            # * one has outgroups and ingroups don't overlap
            merge = False
            if (oa and ob) and not (ia or ib):
                merge = True
            elif (ia or ib) and not (oa or ob):
                merge = True
            elif ((oa and not ob) or (ob and not oa)) and not overlap:
                merge = True

            # print node_id, a, b, "oa=",oa, "ob=", ob, "ia=", ia, "ib=",ib,
            # "ovl=",overlap, "merge=",merge

            if merge:
                nout[node_id] = sum([nout[x] for x in node.succ])
                nin[node_id] = sum([nin[x] for x in node.succ])
                taxa[node_id] = taxa[a].union(taxa[b])
                otus[node_id] = otus[a].union(otus[b])
            else:
                if ia and oa and ib and ob:
                    # write two complete subtrees
                    nout[node_id] = 0
                    nin[node_id] = 0
                    taxa[node_id] = set()
                    otus[node_id] = set()
                    clusters.append(otus[a])
                    clusters.append(otus[b])
                elif ia and oa:
                    # write a, keep b
                    nout[node_id] = nout[b]
                    nin[node_id] = nin[b]
                    taxa[node_id] = taxa[b]
                    otus[node_id] = otus[b]
                    clusters.append(otus[a])
                elif ib and ob:
                    # write b, keep a
                    nout[node_id] = nout[a]
                    nin[node_id] = nin[a]
                    taxa[node_id] = taxa[a]
                    otus[node_id] = otus[a]
                    clusters.append(otus[b])
                elif not(ia and ib and oa and ob):
                    # two empty subtrees merge
                    nout[node_id] = 0
                    nin[node_id] = 0
                    taxa[node_id] = set()
                    otus[node_id] = set()
                else:
                    tree.display()
                    print node_id, ia, ib, oa, ob
                    raise "sanity check failed: unknown case."

    TreeTools.TreeDFS(tree, tree.root,
                      post_function=update_groups)

    # special treatment of root
    if otus[tree.root]:
        if clusters:
            # add to previous cluster if only outgroups or ingroups
            # in root node
            oa = nout[tree.root] > 0
            ia = nin[tree.root] > 0
            if (oa and not ia) or (ia and not oa):
                clusters[-1] = clusters[-1].union(otus[tree.root])
            else:
                clusters.append(otus[tree.root])
        else:
            clusters.append(otus[tree.root])

    return clusters

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printDuplications(duplications, species_nexus, options):

    outfile, is_new = getFile(options, "duplications")
    outfile.write("%s\ttaxa\n" % (Duplication().getHeader()))

    for duplication in duplications:
        outfile.write("%s\t%s\n" % (str(duplication),
                                    ",".join(species_nexus.trees[duplication.mSpeciesTree].get_taxa(duplication.mSpeciesNode))))

    if outfile != options.stdout:
        outfile.close()


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def readDuplications(infile):
    """read duplications from infile from previous analysis.
    """

    duplications = []
    first = True
    for line in infile:
        if first:
            first = False
            continue

        duplication = Duplication()
        duplication.readFromString(line)

        duplications.append(duplication)

    return duplications
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def getReferenceHeight(distance_to_root, gene_tree, input_set,
                       options, extract_species=str,
                       method="median"):
    """return height of tree.

    applies filtering opions from options. Method describes how the
    gene tree is computed from the leaf distances. Possible values are
    median, mean, min, max

    """

    # get maximum height - only genes that are not outgroups
    if options.outgroup_species:
        max_length_set = filter(lambda x: extract_species(
            gene_tree.node(x).data.taxon) not in options.outgroup_species, input_set)
    else:
        max_length_set = input_set

    if len(max_length_set) == 0:
        return None

    if method == "median":
        reference_height = numpy.median(
            map(lambda x: distance_to_root[x], max_length_set))
    elif method == "min":
        reference_height = min(
            map(lambda x: distance_to_root[x], max_length_set))
    elif method == "max":
        reference_height = max(
            map(lambda x: distance_to_root[x], max_length_set))
    elif method == "mean":
        reference_height = numpy.mean(
            map(lambda x: distance_to_root[x], max_length_set))
    else:
        raise "Unknown method %s" % method

    if reference_height == 0.0:
        return None

    return reference_height

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def addCounts(dest, src):
    """add counts in src to dest.
    """
    if type(src) == ListType:
        if src:
            if type(src[0]) not in (DictType, ListType):
                for x in range(len(src)):
                    dest[x] += src[x]
            else:
                for x in range(len(src)):
                    addCounts(dest[x], src[x])

    elif type(src) == DictType:

        for a, b in src.items():
            if type(b) in (FloatType, IntType):
                dest[a] += src[a]
            else:
                addCounts(dest[x], src[x])
    else:
        dest += src

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def appendCounts(dest, src):
    """add counts in src to dest.
    """
    if type(src) == ListType:
        if src:
            if type(src[0]) not in (DictType, ListType):
                dest += src
            else:
                for x in range(len(src)):
                    appendCounts(dest[x], src[x])

    elif type(src) == DictType:

        for x in src.keys():
            if x not in dest:
                if type(src[x]) == DictType:
                    dest[x] = {}
                elif type(src[x]) == ListType:
                    dest[x] = []

            appendCounts(dest[x], src[x])
    else:
        raise "what"
