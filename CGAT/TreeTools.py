"""TreeTools.py - Tools for working with trees
===========================================

This module contains functions to work with gene and/or species trees.

Reference
---------

"""

import Bio
import Bio.Nexus.Nexus
import Bio.Nexus.Trees
import string
import re
from six import StringIO
from CGAT import Tree as Tree


def Newick2Nexus(infile):
    """convert newick formatted tree(s) into a nexus object.

    Multiple trees are separated by a semicolon. Tree names can
    be given by fasta-style separators, i.e., lines starting with
    '>'.

    If the token [&&NHX is found in the tree, it is assumed to be
    output from njtree and support values are added. Support values are
    added in the format taxon:support:branchlength

    Arguments
    ---------
    infile : object
        Input data. Can be a file, a list of lines or a single line.

    Returns
    -------
    nexus : Bio.Nexus.Nexus

    """
    lines = ["#NEXUS\nBegin trees;\n"]

    if isinstance(infile, str):
        tlines = [infile]
    else:
        tlines = [x for x in infile]

    f = []
    id = None
    ntrees = 0

    def __addLine(id, f, lines):

        if len(f) == 0:
            return
        if not id:
            id = "tree%i" % ntrees
        id = re.sub("=", "_", id)

        s = "".join(f)[:-1]
        if s.find("[&&NHX") >= 0:
            # process njtree trees with bootstrap values
            fragments = []
            l = 0
            for x in re.finditer("(:[-0-9.]+\[&&NHX[^\]]*\])", s):
                fragments.append(s[l:x.start()])
                frag = s[x.start():x.end()]
                bl = ":%s" % re.search(":([-0-9.]+)", frag).groups()[0]

                rx = re.search("B=([-0-9.]+)", frag)
                if rx:
                    support = ":%s" % rx.groups()[0]
                else:
                    support = ""

                fragments.append("%s%s" % (support, bl))

                l = x.end()

            fragments.append(s[l:len(s)])

            s = "".join(fragments)
            s = re.sub("\[&&NHX[^\]]*\]", "", s)

        lines.append("tree '%s' = %s;\n" % (id, s))

    for line in tlines:
        line = line.strip()
        if not line:
            continue
        if line[0] == "#":
            continue
        if line[0] == ">":
            id = line[1:]
            continue

        line = re.sub("\s", "", line).strip()
        f.append(line)
        if line[-1] == ";":
            __addLine(id, f, lines)
            f, id = [], None
            ntrees += 1

    # treat special case of trees without trailing semicolon
    __addLine(id, f, lines)
    lines.append("End;")

    # previoulsy, a string was ok, now a string
    # is interpreted as a filename.
    nexus = Bio.Nexus.Nexus.Nexus(StringIO("".join(lines)))

    if len(nexus.trees) == 0:
        raise ValueError("no tree found in file %s" % str(infile))

    # remove starting/ending ' from name
    for tree in nexus.trees:
        tree.name = tree.name[1:-1]

    Tree.updateNexus(nexus)

    return nexus


def Nexus2Newick(nexus, with_branchlengths=True, with_names=False,
                 write_all_taxa=False):
    """convert nexus tree format to newick format.

    Arguments
    ---------
    nexus : Bio.Nexus.Nexus
         The trees to output
    with_branch_lengths : bool
         If True, output branchlengths.
    with_names : bool
         If True, add node names.
    write_all_taxa : bool
         Ouput taxa for internal nodes.

    Returns
    -------
    output : string
         Trees in Newick format.
    """
    lines = []

    for tree in nexus.trees:
        if with_names:
            lines.append(">%s" % tree.name)

        lines.append(Tree2Newick(tree,
                                 with_branch_lengths=with_branchlengths,
                                 write_all_taxa=write_all_taxa))

    return "\n".join(lines)


def Tree2Newick(tree, with_branch_lengths=True, write_all_taxa=False):
    """convert tree to newick format."""
    s = tree.to_string(branchlengths_only=with_branch_lengths,
                       write_all_taxa=write_all_taxa)
    return string.strip(s.split("=")[1])


def Newick2Tree(txt):
    """convert tree to nexus format."""
    return Newick2Nexus(txt).trees[0]


def WriteNexus(nexus, **kwargs):
    """write trees in nexus file format.
    """
    lines = ["#NEXUS\nBegin trees;\n"]
    for t in nexus.trees:
        lines.append("%s\n" % (t.to_string(**kwargs)))
    lines.append("\nEnd;\n")
    return "\n".join(lines)


def GetTaxa(tree):
    """retrieve all taxa of leaves in a tree."""
    return [tree.node(x).get_data().taxon for x in tree.get_terminals()]


def GetTaxonomicNames(tree):
    """get list of taxa."""
    return GetTaxa(tree)


def MapTaxa(tree, map_old2new, remove_unknown=False):
    """update taxa in tree to new taxa.

    Arguments
    ---------
    tree : Tree
        The tree to update.
    map_old2new :dict
        Dictionary mapping old taxa to new taxa.
    remove_unknown : bool
        If true, taxa not in `map_old2new` will be removed.
    """

    unknown = []
    for n, node in list(tree.chain.items()):
        if node.data.taxon:
            try:
                node.data.taxon = map_old2new[node.data.taxon]
            except KeyError:
                unknown.append(node.data.taxon)

    if remove_unknown:
        for taxon in unknown:
            tree.prune(taxon)


def Branchlength2Support(tree):
    """copy values stored as branchlength to into support

    The branchlength property is not changed.

    This step is necessary when support has been stored as branchlength
    (e.g. paup), and has thus been read in as branchlength.
    """

    for n in list(tree.chain.keys()):
        tree.node(n).data.support = tree.node(n).data.branchlength


def Species2Genes(nexus, map_species2genes):
    """convert a species tree to a gene tree.

    Arguments
    ---------
    nexus : Bio.Nexus.Nexus
         The trees to work on
    map_species2genes : dict
         Dictionary mapping species names to gene names
    """

    for tree in nexus.trees:
        for nx in tree.get_terminals():
            t1 = tree.node(nx).get_data().taxon
            if t1 in map_species2genes:
                for g in map_species2genes[t1]:
                    d = Nexus.NodeData(taxon=g)
                    new_node = Node(d)
                    tree.add(new_node, nx)
                    tree.node(nx).get_data().taxon = None


def Genes2Species(nexus, map_gene2species):
    """convert a gene tree into a species tree.

    Arguments
    ---------
    nexus : Bio.Nexus.Nexus
         The trees to work on
    map_gene2species : dict
         Dictionary mapping gene names to species names

    """

    for tree in nexus.trees:
        for nx in tree.get_terminals():
            t1 = tree.node(nx).get_data().taxon
            if t1 in map_gene2species:
                tree.node(nx).get_data().taxon = map_gene2species[t1]


def BuildMapSpecies2Genes(genes, pattern_species="^([^|]+)[|]"):
    """build a map of species to genes

    This method assumes that gene names contain the species name and
    it can be extracted via a regular expression.

    Arguments
    ---------
    genes : list
         List of genes
    pattern_species : string
         Regular expression to extract species name from gene name.

    Returns
    -------
    map_species2genes : dict
         Mapping between species to one or more genes
    map_gene2species : dict
         Mapping between a gene to the species
    """
    rx = re.compile(pattern_species)
    map_species2genes = {}
    map_gene2species = {}
    for gene in genes:
        species = rx.search(gene).groups()[0]
        if species not in map_species2genes:
            map_species2genes[species] = []
        map_species2genes[species].append(gene)
        map_gene2species[gene] = species

    return map_species2genes, map_gene2species


def GetMonophyleticPairs(tree):
    """build list of monophyletic pairs in tree.
    """

    leaves = tree.get_terminals()

    pairs = []

    for z in range(len(leaves) / 2):

        for x in range(0, len(leaves) - 1):
            rx = leaves[x]
            tx = tree.node(rx).get_data().taxon
            for y in range(x + 1, len(leaves)):
                ry = leaves[y]
                ty = tree.node(ry).get_data().taxon
                if tree.is_monophyletic((tx, ty)) != -1:
                    pairs.append((rx, ry, tx, ty))

    return pairs


def GetTaxaForSpecies(tree, species, pattern_species="^([^|]+)[|]"):
    """get all taxa of a given species.

    This method assumes that node labels contain the species name and
    it can be extracted via a regular expression.

    Arguments
    ---------
    genes : list
         List of genes
    pattern_species : string
         Regular expression to extract species name from gene name.

    Returns
    -------
    taxa : list
         List of taxa from this species.

    """
    rx = re.compile(pattern_species)

    taxa = []
    for l in GetTaxa(tree):
        g = rx.search(l).groups()[0]
        if g == species:
            taxa.append(l)

    return taxa


def IsMonophyleticForSpecies(tree,
                             species,
                             pattern_species="^([^|]+)[|]"):
    """check if a tree is monophyletic for a species.

    This method assumes that node labels contain the species name and
    it can be extracted via a regular expression.

    Arguments
    ---------
    tree : :class:`Tree`
         Tree to analyse
    species : string
         Species to check
    pattern_species : string
         Regular expression to extract species name from gene name.

    Returns
    -------
    bool
    """

    taxa = GetTaxaForSpecies(tree, species, pattern_species)
    tree.root_with_outgroup(taxa)
    return tree.is_monophyletic(taxa) != -1


def IsMonophyleticForTaxa(tree,
                          taxa,
                          support=None):
    """check if a tree is monophyletic for a list of taxa.

    Arguments
    ---------
    tree : :class:`Tree`
         Tree to analyse
    taxa : list
         List of taxa
    support : float
         Minimum bootstrap support

    Returns
    -------
    bool

    """
    tree.root_with_outgroup(taxa)

    if support:
        n = tree.is_monophyletic(taxa)
        if n == -1:
            return False
        return tree.node(tree.node(tree.root).succ[0]).data.support >= support
    else:
        return tree.is_monophyletic(taxa) != -1


def GetLeaves(tree, node):
    """Return leaves in tree below node.
    """
    tree.root_with_outgroup((node,))
    return tree.get_taxa(node)


def IsSingleSpecies(tree, node, pattern_species="^([^|]+)[|]"):
    """True if taxa below node contain the same species."""
    rx = re.compile(pattern_species)

    species = {}
    leaves = GetLeaves(tree, node)
    for l in leaves:
        species[rx.search(l).groups()[0]] = 1

    return len(species) == 1


def CountDuplications(tree, species,
                      pattern_species="^([^|]+)[|]"):
    """Count the number duplications for a given species.

    Do not check for monophyly versus species.
    """
    result = []
    taxa = GetTaxaForSpecies(tree, species, pattern_species)
    ids = list(map(tree.search_taxon, taxa))
    for x in range(len(ids) - 1):
        for y in range(x + 1, len(ids)):
            t = tree.common_ancestor(ids[x], ids[y])
            is_single_species = IsSingleSpecies(tree, t, pattern_species)
            tree.root_with_outgroup((taxa[x], taxa[y]))
            is_monophyletic = tree.is_monophyletic((taxa[x], taxa[y])) != -1
            result.append(
                (taxa[x], taxa[y], is_single_species,
                 is_monophyletic, tree.distance(ids[x], ids[y])))

    return result


def Transcript2GeneTree(tree,
                        map_transcript2gene,
                        map_gene2transcripts):
    """convert a transcript tree into a gene tree.

    supply a map for mapping transcripts to genes.

    The procedure for converting a transcript tree into a gene tree:

    If there are two genes, and they are monophyletic, no matter how many
    transcripts, the order is as follows:

       1 Merge all nodes into two, one for each gene.

       2 The distance between the genes is the minimum distance observed between
         two transcripts from different genes. Half of this will be set as the
         branch length from the gene leaves.

    If this is not possible for a set of genes, the procedure will fail and not
    return a gene tree.
    """
    raise NotImplementedError()
    MapTaxa(tree, map_transcript2gene)

    # get all leaves and sort by taxon
    ids = tree.get_terminals()

    # sort identities by taxa
    ids.sort(
        lambda x, y: cmp(tree.node(x).get_data().taxon, tree.node(y).get_data().taxon))
    taxa = [tree.node(x).get_data().taxon for x in ids]

    print(ids)
    print(taxa)

    for x in range(len(taxa) - 1):
        for y in range(x + 1, len(taxa)):
            print(ids[x], ids[y], taxa[x], taxa[y])


def MapTerminalTaxa(tree, mapping):
    """map taxa in leaves in all trees."""

    for nx in tree.get_terminals():
        t1 = tree.node(nx).get_data().taxon
        if t1 in mapping:
            tree.node(nx).get_data().taxon = mapping[t1]


def GetCommonAncestor(tree, taxa):
    """retrieve common ancestor for a list of taxa.

    Reroot tree. Check if it is monopyletic. If it is, return root,
    otherwise, return -1.
    """

    tree.root_with_outgroup(taxa)
    if tree.is_monophyletic(taxa) == -1:
        return -1
    else:
        x = tree.search_taxon(taxa[0])
        trace = tree.trace(tree.root, x)

        if trace[0] == tree.root:
            return tree.root
        else:
            return trace[0]


def Nop(x):
    return True


def TreeDFS(tree, node_id,
            pre_function=Nop,
            descend_condition=Nop,
            post_function=Nop):
    """BFS tree tree traversal starting at node_id.

    Apply functions pre_function at first and
    post_function at last visit of a node.
    """
    pre_function(node_id)
    for n in tree.node(node_id).succ:
        if descend_condition(node_id):
            TreeDFS(tree, n, pre_function, descend_condition, post_function)
    post_function(node_id)


def GetMaxIndex(tree):
    """get maximum node number."""
    return tree.id


def GetNumChildren(tree):

    nnodes = GetMaxIndex(tree) + 1

    counts = [0] * nnodes

    def count(node_id):
        s = tree.node(node_id).succ
        if s == []:
            counts[node_id] = 1
        else:
            for n in s:
                counts[node_id] += counts[n]

    TreeDFS(tree, tree.root, post_function=count)
    return counts


def GetBranchLengths(tree):
    """return an array with minimum and maximum branch length."""

    nnodes = GetMaxIndex(tree) + 1

    max_sums = [0] * nnodes
    min_sums = [0] * nnodes

    def count(node_id):
        s = tree.node(node_id).succ
        if s != []:
            min_vals = [
                min_sums[n] + tree.node(n).data.branchlength for n in s]
            max_vals = [
                max_sums[n] + tree.node(n).data.branchlength for n in s]

            min_sums[node_id] = min(min_vals)
            max_sums[node_id] = max(max_vals)

    TreeDFS(tree, tree.root, post_function=count)

    return min_sums, max_sums


def Reroot(tree, taxa):
    """reroot tree with taxa - the list of
    taxa does not need to be monophyletic.
    """

    nnodes = GetMaxIndex(tree) + 1

    within_taxa = [False] * nnodes

    def update_taxa(node_id):
        n = tree.node(node_id)

        s = n.succ
        if s == []:
            within_taxa[node_id] = n.data.taxon in taxa
        else:
            for ss in s:
                within_taxa[node_id] |= within_taxa[ss]

    TreeDFS(tree, tree.root,
            post_function=update_taxa)

    # go down the tree - get largest branch
    # in each subtree left/right from root
    # that contain all matches to taxa
    # Check if root spans taxa
    all_true = True
    for n in tree.node(tree.root).succ:
        all_true &= within_taxa[n]

    # if root spans taxa get largest subtrees with no matches to taxa
    if all_true:
        def has_taxa(node_id):
            if not within_taxa[node_id]:
                extra_subtree[node_id] = True
                return False
            else:
                return True
    else:
        # if root does not span taxa, get smallest subtree including
        # all taxa
        def has_taxa(node_id):
            """check if all children contain taxa."""
            if not within_taxa[node_id]:
                return False
            for n in tree.node(node_id).succ:
                if not within_taxa[n]:
                    return True
            extra_subtree[node_id] = True
            return False

    # do the search
    extra_subtree = [False] * nnodes
    TreeDFS(tree, tree.root,
            descend_condition=has_taxa)

    nodes = [x for x in range(nnodes) if extra_subtree[x]]

    if len(nodes) == 0:
        # no rerooting, if all or no taxa within tree
        return tree.root
    elif len(nodes) > 1:
        # if more than two nodes (i.e, if all_true is True)
        nchildren = GetNumChildren(tree)
        nodes.sort(lambda x, y: cmp(nchildren[x], nchildren[y]))
        nodes.reverse()

    subtree_node = nodes[0]

    # todo: write rerooting function from node
    taxa = tree.get_taxa(subtree_node)
    result = tree.root_with_outgroup(taxa)

    if result == -1:
        raise TreeError('error while rerooting tree')

    # return node_id with all taxa
    if all_true:
        for n in tree.node(tree.root).succ:
            if n != subtree_node:
                return n
    else:
        return subtree_node


def GetSubsets(tree, node=None, with_decoration=True):
    """return subsets below a certain node including
    their height (distance from leaves) and branchlength
    """
    if node is None:
        node = tree.root

    if with_decoration:
        if tree.node(node).succ == []:
            return [((tree.node(node).data.taxon,), 0, tree.node(node).data.branchlength,)]
        else:
            # get all subtrees
            children, height, branchlength = [], 0, tree.node(
                node).data.branchlength
            subtrees = [GetSubsets(tree, n) for n in tree.node(node).succ]

            ss = []
            for s in subtrees:
                children += s[-1][0]
                height += s[-1][1] + s[-1][2]
                ss += s

            height /= len(tree.node(node).succ)

            return ss + [(children, height, branchlength)]
    else:
        if tree.node(node).succ == []:
            return [[tree.node(node).data.taxon, ]]
        else:
            # get all subtrees
            subtrees = [GetSubsets(tree, n, with_decoration=False)
                        for n in tree.node(node).succ]

            ss = []
            children = []
            for s in subtrees:
                children += s[-1]
                ss += s

            return ss + [children]


def CountBranchPoints(tree, taxa):
    """count the number branch points together with their
    distances for a given list of taxa.

    return a list of branch points
    """
    tree.root_with_outgroup(taxa)

    if tree.is_monophyletic(taxa) == -1:
        return None

    parent = GetCommonAncestor(tree, taxa)

    # retrieve all subsets with their branchlengths.
    return GetSubsets(tree, parent)


def IsCompatible(tree1, tree2):
    """check if two trees are compatible.

    note: this will delete support information.
    """

    if len(tree1.get_terminals()) != len(tree2.get_terminals()):
        return False, "leaves"

    l1 = GetTaxonomicNames(tree1)
    l2 = GetTaxonomicNames(tree2)

    for n in list(tree2.chain.keys()):
        tree2.node(n).data.support = 1.0
    for n in list(tree1.chain.keys()):
        tree1.node(n).data.support = 1.0

    for l in l1:
        if l not in l2:
            return False, "taxa"

    if tree1.is_compatible(tree2, 0):
        return False, "topology"
    else:
        return True, ""


def Tree2Graph(tree):
    """return tree as a list of edges in a graph."""
    links = []
    for node_id1, node1 in list(tree.chain.items()):
        if node1.prev is not None:
            links.append((node_id1, node1.prev, node1.get_data().branchlength))
    return links


def Graph2Tree(links, label_ancestral_nodes=False):
    """build tree from list of nodes.

    Assumption is that links always point from parent to child.
    """
    tree = Tree.Tree()

    # map of names to nodes in tree
    map_node2id = {links[0][0]: 0}
    map_id2node = [links[0][0]]

    for parent, child, branchlength in links:
        if parent not in map_node2id:
            p = len(tree.chain)
            tree.chain[p] = Bio.Nexus.Nodes.Node(Bio.Nexus.Trees.NodeData())
            map_node2id[parent] = p
            map_id2node.append(parent)
        else:
            p = map_node2id[parent]

        if child not in map_node2id:
            c = len(tree.chain)
            tree.chain[c] = Bio.Nexus.Nodes.Node(Bio.Nexus.Trees.NodeData())
            map_node2id[child] = c
            map_id2node.append(child)
        else:
            c = map_node2id[parent]

        tree.chain[p].succ.append(c)
        tree.chain[c].prev = p
        tree.chain[c].data.branchlength = branchlength

    # set taxon names for children and find root
    for i, n in list(tree.chain.items()):
        if n.prev == []:
            tree.root = i

        if n.succ == [] or label_ancestral_nodes:
            n.data.taxon = map_id2node[i]

    # set pointer to last id
    tree.id = len(list(tree.chain.items())) - 1

    return tree


def GetAllNodes(tree):
    """return all nodes in the tree."""
    return list(tree.chain.keys())


def GetDistancesBetweenTaxa(tree, taxa1, taxa2):
    """get average branchlength between taxa1 and taxa2."""

    # get sets with terminal nodes in taxa
    a, b = [], []
    for x in tree.get_terminals():
        t = tree.node(x).data.taxon
        if t in taxa1:
            a.append((t, x))
        elif t in taxa2:
            b.append((t, x))

    distances = []
    for ta, aa in a:
        for tb, bb in b:
            distances.append((ta, tb, tree.distance(aa, bb)))

    return distances


def PruneTerminal(tree, taxon):
    """Prunes a terminal taxon from the tree.

    id_of_previous_node = prune(tree,taxon)
    If taxon is from a bifurcation, the connecting node will be collapsed
    and its branchlength added to remaining terminal node. This might be no
    longer a meaningful value.

    direct copy of Nexus.Trees.py - don't know why have a separate method,
    maybe there was a bug in Nexus.Trees.
    """
    return tree.prune(taxon)

    id = tree.search_taxon(taxon)
    if id is None:
        raise ValueError('Taxon not found: %s' % taxon)
    elif id not in tree.get_terminals():
        raise ValueError('Not a terminal taxon: %s' % taxon)
    else:
        prev = tree.unlink(id)
        tree.kill(id)
        if not prev == tree.root and len(tree.node(prev).succ) == 1:
            succ = tree.node(prev).get_succ()
            tree.collapse(prev)

        return prev


def add_children(old_tree, new_tree, old_id, new_id):

    for n in old_tree.node(old_id).succ:
        nid = new_tree.add(Bio.Nexus.Nodes.Node(old_tree.node(n).data), new_id)
        add_children(old_tree, new_tree, n, nid)


def GetSubtree(tree, node_id):
    """return a copy of tree from node_id downwards."""

    subtree = Tree.Tree(weight=tree.weight,
                        rooted=tree.rooted,
                        name=tree.name)

    # automatically adds a root, so substitute it
    n = Bio.Nexus.Nodes.Node(tree.node(node_id).data)
    n.id = subtree.root
    subtree.chain[subtree.root] = n

    add_children(tree, subtree, node_id, subtree.root)
    return subtree


def Unroot(tree):
    """unroot tree."""

    # collapse a child of the root, that has
    # at least two children and where both
    # children are not terminals.
    root_node = tree.node(tree.root)
    root_id = tree.root

    # check if root has a single child - if so, remove this root
    if len(root_node.succ) == 1:
        root_id = root_node.succ[0]
        tree.kill(tree.root)
        tree.root = root_id
        root_node = tree.node(tree.root)

    if len(tree.node(root_node.succ[0]).succ) == 0 and len(tree.node(root_node.succ[1]).succ) == 0:
        return

    # calculate branch length along branch that has
    # been split by root. This value needs to be assigned
    # to the remaining childs branchlength.
    n = sum([tree.node(x).data.branchlength for x in root_node.succ])

    if len(tree.node(root_node.succ[0]).succ) > 1:
        x = root_node.succ[0]
        y = root_node.succ[1]
    else:
        x = root_node.succ[1]
        y = root_node.succ[0]

    tree.collapse(x)
    tree.node(y).data.branchlength = n
    tree.rooted = False


def GetSize(tree):
    """return the length of the tree. This is the maximum node_id + 1.

    This quantity is useful for tree traversal while updating
    a container.
    """
    return max(tree.chain.keys()) + 1


def PruneTree(tree, taxa, keep_distance_to_root=False):
    """prune tree: keep only those taxa in list.
    """

    for nx in tree.get_terminals():
        taxon = tree.node(nx).get_data().taxon
        if taxon not in taxa:
            tree.prune(taxon)

    r = tree.root
    rn = tree.node(tree.root)

    # if one complete side of the root has been removed,
    # collapse it.
    if len(rn.succ) == 1:
        # tree.collapse on root does not work
        # tree.collapse( tree.root )
        s = rn.succ[0]
        sn = tree.node(s)
        sn.prev = None
        if not keep_distance_to_root:
            sn.data.branchlength = 0.0
        tree.root = s
        tree.kill(r)


def GetNodeMap(tree1, tree2):
    """map nodes between tree1 and tree2.
    """

    if not tree1.is_identical(tree2):
        raise ValueError("trees are not the same")

    map_a2b = [0] * len(list(tree1.chain.keys()))

    for i, n in list(tree1.chain.items()):
        t = tree1.get_taxa(i)
        o = tree2.is_monophyletic(t)
        if o != -1:
            map_a2b[i] = o
        else:
            raise ValueError("trees are not congruent.")

    return map_a2b


class NodeType:
    mType = "Generic"
    mDescription = "Generic"

    def __init__(self, node1, node2):
        self.mSpeciesNode = node1
        self.mGeneNode = node2

    def __str__(self):
        return "\t".join((self.mType,
                          str(self.mSpeciesNode),
                          str(self.mGeneNode)))


class NodeTypeSpeciation(NodeType):

    mType = "Speciation"
    mDescription = "Speciation event"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeSpeciationDeletion(NodeType):

    mType = "SpeciationDeletion"
    mDescription = "Speciation event, but in one sub-branch, "
    "deletions occured."

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeDuplication(NodeType):
    mType = "Duplication"
    mDescription = "Duplication event"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeOutparalogs(NodeType):
    mType = "Outparalogs"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeDuplicationDeletion(NodeType):
    mType = "DuplicationDeletion"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeDuplicationLineage(NodeType):
    mType = "DuplicationLineage"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeDuplicationInconsistency(NodeType):
    mType = "DuplicationInconsistency"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeTranscripts(NodeType):
    mType = "Transcripts"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeInconsistency(NodeType):
    mType = "Inconsistency"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeInconsistentTranscripts(NodeType):
    mType = "InconsistentTranscripts"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeMasked(NodeType):
    mType = "Masked"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


class NodeTypeLeaf(NodeType):
    mType = "Leaf"

    def __init__(self, *args, **kwargs):
        NodeType.__init__(self, *args, **kwargs)


def ReconciliateByRio(gene_tree, species_tree,
                      extract_species,
                      extract_gene=None,
                      outgroup_species=None,
                      min_branch_length=0.0):
    """
    Gene tree G and species tree S

    If outgroup_species is given: trees will be cut of
    as soon as one of the outgroup species is part of a subtree.
    The corresponding node type will be out-paralog. Out-paralog
    relationship is cast upwards.

    Input trees are rooted and binary.

    Output: gene tree with duplication/speciation assigned to each node.

    Initialization:

        Number nodes in S in pre-order traversal (root = 1), such
        that child nodes are always larger than parent nodes.

        For each external node g of G, set M(g) to the number of the
        external node in S with the matching species name.

    Recursion:

        Visit each internal node g of G in post-order traversal, (i.e.
        from leaves to root)::

          set a = M(g1) # g1 = first child of current node g
          set b = M(g2) # g2 = second child of current node g

          while a != b:
              if a > b:
                    set a = parent of node a in species tree
              else:
                    set b = parent of node b in species tree
          set M(g) = a

          if M(g) == M(g1) or M(g) == M(g2):
              g is duplication
          else:
              g is speciation

    The algorithm returns an array for each node with its type. 

    If extract_gene is given, the algorithm will label transcription nodes
    for alternative transcripts (duplications involving the same gene).

    The algorithm has been extended to accomodate the following test cases:

    Alternative transcripts
        Alternative transcripts that span genes from other species are permitted,
        if at most one gene of the other species is involved.

        To avoid over-counting of speciation events, the one subtree with the
        least species is masked.

    If the branch length of a node in the gene tree is shorter than min_branch_length,
    the resultant node is masked, because the topology might be dodgy.
    """

    ########################################################################
    # Initialization
    nnodes_genetree = max(gene_tree.chain.keys()) + 1
    nnodes_speciestree = max(species_tree.chain.keys()) + 1
    nspecies = len(species_tree.get_taxa())

    # vector with node numbers
    N = [0] * (nnodes_speciestree + 1)
    map_N2Parent = [0] * (nnodes_speciestree + 2)
    map_N2node_id = [0] * (nnodes_speciestree + 2)
    map_species2N = {}
    # Mapping function
    M = [0] * nnodes_genetree
    # result: speciation array
    node_types = [None] * nnodes_genetree

    m = N[0]

    # relabel the species tree externally
    counter = [1]

    def init(node_id):
        node = species_tree.node(node_id)
        if node.succ == []:
            map_species2N[node.data.taxon] = counter[0]
        map_N2node_id[counter[0]] = node_id
        N[node_id] = counter[0]
        if node.prev:
            map_N2Parent[counter[0]] = N[node.prev]
        counter[0] += 1

    TreeDFS(species_tree, species_tree.root,
            pre_function=init)

    if outgroup_species:
        outgroups = set(outgroup_species)
    else:
        outgroups = None

    for x in gene_tree.get_terminals():
        M[x] = map_species2N[extract_species(gene_tree.node(x).data.taxon)]

    ########################################################################
    ########################################################################
    # Recursion for masking a subtree
    def mask_subtree(node_id):
        node = gene_tree.node(node_id)
        if node.succ == []:
            return

        node_types[node_id] = NodeTypeMasked(node_types[node_id].mSpeciesNode,
                                             node_types[node_id].mGeneNode)

    ########################################################################
    ########################################################################
    # Recursion for updating assignments
    def update(g):

        ng = gene_tree.node(g)
        if ng.succ == []:
            node_types[g] = NodeTypeLeaf(map_N2node_id[M[g]], g)
            return

        if len(ng.succ) != 2:
            gene_tree.display()
            raise ValueError("warning: not a binary tree.")

        g1, g2 = ng.succ
        a = M[g1]
        b = M[g2]

        while a != b:
            if a > b:
                a = map_N2Parent[a]
            else:
                b = map_N2Parent[b]

        M[g] = a

        taxa1 = gene_tree.get_taxa(g1)
        taxa2 = gene_tree.get_taxa(g2)
        species1 = set(map(extract_species, taxa1))
        species2 = set(map(extract_species, taxa2))

        if min_branch_length:
            # mask nodes with short branches on this node or to both children
            min_bl = min(ng.data.branchlength, min(
                [gene_tree.node(x).data.branchlength for x in ng.succ]))
            if min_bl < min_branch_length:
                node_types[g] = NodeTypeMasked(map_N2node_id[M[g]], g)
                return

        if M[g] == 1:
            # check if species sets of subtrees are overlapping
            # If they are, then we have outparalogs, otherwise, it is the final
            # speciation event
            if species1.intersection(species2):
                node_types[g] = NodeTypeOutparalogs(map_N2node_id[M[g]], g)
            else:
                node_types[g] = NodeTypeSpeciation(map_N2node_id[M[g]], g)
        elif outgroups and (species1.intersection(outgroups) or species2.intersection(outgroups)):
            node_types[g] = NodeTypeOutparalogs(map_N2node_id[M[g]], g)
        elif M[g] == M[g1] or M[g] == M[g2]:
            # additional check: check if species sets of subtrees are actually
            # overlapping
            if species1.intersection(species2):

                # check for inconsistent transcripts
                if extract_gene:
                    try:
                        genes1 = set(
                            zip(list(map(extract_species, taxa1)), list(map(extract_gene, taxa1))))
                    except AttributeError:
                        raise AttributeError(
                            "could not parse %s" % (",".join(taxa1)))
                    try:
                        genes2 = set(
                            zip(list(map(extract_species, taxa2)), list(map(extract_gene, taxa2))))
                    except AttributeError:
                        raise AttributeError(
                            "could not parse %s" % (",".join(taxa2)))

                    if genes1.intersection(genes2):
                        if len(genes1) == len(genes2) and len(genes1) == 1:
                            # transcripts of the same gene are joined
                            node_types[g] = NodeTypeTranscripts(
                                map_N2node_id[M[g]], g)
                        else:
                            # for all transcripts that come from the same genome - check if
                            # they are all from the same gene(s). If they are, we have an alternative
                            # transcripts node.
                            is_ok = True
                            for species in species1.intersection(species2):
                                sg1 = set(
                                    [x for x in genes1 if x[0] == species])
                                sg2 = set(
                                    [x for x in genes2 if x[0] == species])
                                if sg1.intersection(sg2) != sg1.union(sg2):
                                    is_ok = False
                            if is_ok:
                                # all transcripts are from the same gene
                                # mark node as alternative transcript one
                                node_types[g] = NodeTypeTranscripts(
                                    map_N2node_id[M[g]], g)
                                # mask subtree with fewer species to avoid over-counting
                                # of speciation events
                                if len(species1) > len(species2):
                                    TreeDFS(
                                        gene_tree, g2, post_function=mask_subtree)
                                else:
                                    TreeDFS(
                                        gene_tree, g1, post_function=mask_subtree)
                            else:
                                # transcripts of the same gene are joined, but
                                # they are not monophyletic
                                node_types[g] = NodeTypeInconsistentTranscripts(
                                    map_N2node_id[M[g]], g)

                if node_types[g] is None:
                    if len(species1) == len(species2) and len(species1) == 1:
                        # lineage specific duplication: only one species is
                        # involved
                        node_types[g] = NodeTypeDuplicationLineage(
                            map_N2node_id[M[g]], g)
                    elif len(species1) == len(species2) and not species1.difference(species2):
                        # clean duplication: species sets do completely overlap
                        node_types[g] = NodeTypeDuplication(
                            map_N2node_id[M[g]], g)
                    elif len(species1.difference(species2)) == 0 or len(species2.difference(species1)) == 0:
                        # duplication, but with only deletions in one branch
                        node_types[g] = NodeTypeDuplicationDeletion(
                            map_N2node_id[M[g]], g)
                    else:
                        node_types[g] = NodeTypeDuplicationInconsistency(
                            map_N2node_id[M[g]], g)
            # check for speciation events with a deletion
            ##
            else:
                # check for alternative transcripts.
                if ng.data.branchlength < min_branch_length:
                    node_types[g] = NodeTypeMasked(map_N2node_id[M[g]], g)
                else:
                    node_types[g] = NodeTypeInconsistency(
                        map_N2node_id[M[g]], g)
        else:
            # get species expected under species tree
            expected_species = set(species_tree.get_taxa(map_N2node_id[M[g]]))
            observed_species = species1.union(species2)
            if len(observed_species.intersection(expected_species)) < len(expected_species):
                node_types[g] = NodeTypeSpeciationDeletion(
                    map_N2node_id[M[g]], g)
            else:
                node_types[g] = NodeTypeSpeciation(map_N2node_id[M[g]], g)

    TreeDFS(gene_tree, gene_tree.root,
            post_function=update)

    return node_types


def CountDuplications(gene_tree, species_tree, node_types,
                      extract_species,
                      extract_gene=None):
    """count duplications.

    given are gene and species tree and node types (duplication/speciation)

    extract_species gives the species for an OTU in the gene tree

    Extract_gene gives the gene for an OTU in the gene tree. If not given,
    all transcripts are counted as unique.
    """

    ########################################################################
    # Get additional data for branch lengths
    min_branch_lengths, max_branch_lengths = GetBranchLengths(gene_tree)

    for x in range(len(node_types)):

        node_type = node_types[x]
        if node_type:
            if node_type.mType not in ("Speciation"):
                print("\t".join(map(str, (node_type,
                                          (min_branch_lengths[node_type.mGeneNode] + max_branch_lengths[node_type.mGeneNode]) / 2))))
                for s in species_tree.node(node_type.mSpeciesNode).succ:
                    print("\t" * 5 + ",".join(species_tree.get_taxa(s)))
                for s in gene_tree.node(node_type.mGeneNode).succ:
                    print("\t" * 5 + ",".join(gene_tree.get_taxa(s)))


def GetParentNodeWhereTrue(node_id, tree, stop_function):
    """walk up in gene tree and stop where stop_function is true.

    The walk finishes at the root.

    returns tuple of node and distance.
    """

    node = tree.node(node_id)
    distance = 0
    while node.prev is not None:

        if stop_function(node_id):
            return node_id, distance

        distance += node.data.branchlength

        node_id = node.prev
        node = tree.node(node_id)

    return node_id, distance


def GetChildNodesWhereTrue(node_id, tree, stop_function):
    """walk down in tree and stop where stop_function is true

    The walk finishes at the leaves.

    returns a list of tuples of nodes and distance.
    """

    result = []

    def __getChildNodes(node_id, distance):

        node = tree.node(node_id)

        distance += node.data.branchlength

        if not node.succ:
            result.append((node_id, distance))
        elif stop_function(node_id):
            result.append((node_id, distance))
        else:
            for x in node.succ:
                __getChildNodes(x, distance)

    node = tree.node(node_id)
    __getChildNodes(node_id, -node.data.branchlength)

    return result


def GetDistanceToRoot(tree):
    """return list with distance to root for each node."""
    #########################################################################
    #########################################################################
    #########################################################################
    # compute distance to root for each node
    #########################################################################
    distance_to_root = [0] * GetSize(tree)

    def record_distance(node_id):
        node = tree.node(node_id)
        if node.prev:
            distance_to_root[
                node_id] += distance_to_root[node.prev] + node.data.branchlength
        else:
            distance_to_root[node_id] = node.data.branchlength

    TreeDFS(tree, tree.root, pre_function=record_distance)

    return distance_to_root


def traverseGraph(graph, start, block=[]):
    """traverse graph, go not passed nodes in block.
    """

    to_visit = [start, ]
    visited = {}

    while to_visit:
        v = to_visit[0]
        del to_visit[0]
        visited[v] = 1
        for n in graph[v]:
            if n not in visited and n not in block:
                to_visit.append(n)

    return visited


def getPattern(tree, nodes, map_taxon2position):

    pattern = ["0"] * len(map_taxon2position)
    for n in nodes:
        t = tree.node(n).get_data().taxon
        if t is not None:
            pattern[map_taxon2position[t]] = "1"
    return "".join(pattern)


def convertTree2Graph(tree):
    """convert tree to a graph."""

    graph = {}
    edges = []
    for i, n in list(tree.chain.items()):
        if i not in graph:
            graph[i] = []
        for nn in n.succ:
            if nn not in graph:
                graph[nn] = []
            graph[nn].append(i)
            graph[i].append(nn)
            edges.append((i, nn))

    return graph, edges


def calculatePatternsFromTree(tree, sort_order):
    """calculate patterns from a tree."""

    notus = len(sort_order)

    map_taxon2position = {}
    for x in range(notus):
        map_taxon2position[sort_order[x]] = x

    graph, edges = convertTree2Graph(tree)
    patterns = []
    for a, b in edges:
        result = traverseGraph(graph, a, [b, ])
        patterns.append(getPattern(tree, list(
            result.keys()), map_taxon2position))
        result = traverseGraph(graph, b, [a, ])
        patterns.append(getPattern(tree, list(
            result.keys()), map_taxon2position))

    patterns.append("1" * notus)
    return patterns
