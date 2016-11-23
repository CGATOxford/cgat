"""
Tree.py - A phylogenetic tree
=============================

The :class:`Tree` is derived from the class from Bio.Nexus.Trees.Tree
adding some additional functionality.

Reference
---------

"""
import math

import Bio.Nexus
import Bio.Nexus.Nodes
import Bio.Nexus.Trees


def updateNexus(nexus):
    """change trees in a nexus object (see Biopython_) to :class:`Tree`.
    """
    for x in range(len(nexus.trees)):
        # remove the trailing ";"
        t = nexus.trees[x]

        # in Biopython, to_string includes the name separated by =, while
        # my own class does not.
        tree = t.to_string(branchlengths_only=False, plain=False)
        a = tree.find("=")
        if a >= 0:
            tree = tree[a + 1:]

        nexus.trees[x] = Tree(tree=tree.strip(),
                              weight=t.weight,
                              rooted=t.rooted,
                              name=t.name,
                              data=t.dataclass,
                              max_support=t.max_support)


def Nop(x):
    """empty function for tree traversal"""
    return True


class Tree(Bio.Nexus.Trees.Tree):
    """A phylogenetic tree.

    This class represents a tree using a chain of nodes with on
    predecessor (=ancestor) and multiple successors (=subclades).
    """

    def __init__(self, *args, **kwargs):
        Bio.Nexus.Trees.Tree.__init__(self, *args, **kwargs)

    def __len__(self):
        """returns the number of nodes in the tree."""
        return len(list(self.chain.keys()))

    def root_at_node(self, node, distance=0):
        """root tree at node.

        Arguments
        ---------
        node :
            New root
        distance : float
            Distance of node to new root.

        This is a subset of the code taken from root_with_outgroup.
        """

        def _connect_subtree(parent, child):
            """Hook subtree starting with node child to parent."""
            for i, branch in enumerate(self.unrooted):
                if parent in branch[:2] and child in branch[:2]:
                    branch = self.unrooted.pop(i)
                    break
            else:
                raise ValueError('Unable to connect nodes for rooting: '
                                 'nodes %d and %d are not connected' %
                                 (parent, child))
            self.link(parent, child)
            self.node(child).data.branchlength = branch[2]
            self.node(child).data.support = branch[3]
            # now check if there are more branches connected to the child, and
            # if so, connect them
            child_branches = [b for b in self.unrooted if child in b[:2]]
            for b in child_branches:
                if child == b[0]:
                    succ = b[1]
                else:
                    succ = b[0]
                _connect_subtree(child, succ)

        # code for getting the outgroup node has been removed. Instead
        # outgroup_node is set to node.
        outgroup_node = node

        if outgroup_node == self.root:
            return self.root

        # find ingroup node on a rooted tree
        ingroup_node = self.node(outgroup_node).prev

        # if tree is already rooted with outgroup on a bifurcating root,
        if ingroup_node == self.root and len(self.node(self.root).succ) == 2:
            s = self.node(ingroup_node).succ
            if s[0] == outgroup_node:
                ingroup_node = s[1]
            else:
                ingroup_node = s[0]

        # move to unrooted structure
        self.unroot()

        # now we find the branch that connects outgroup and ingroup
        # print self.node(outgroup_node).prev
        for i, b in enumerate(self.unrooted):
            if outgroup_node in b[:2] and ingroup_node in b[:2]:
                root_branch = self.unrooted.pop(i)
                break
        else:
            raise ValueError('Unrooted and rooted Tree do not match')

        # now we destroy the old tree structure, but keep node data. Nodes will
        # be reconnected according to new outgroup
        for n in self.all_ids():
            self.node(n).prev = None
            self.node(n).succ = []

        # now we just add both subtrees (outgroup and ingroup) branch for
        # branch
        root = Bio.Nexus.Nodes.Node(
            data=Bio.Nexus.Trees.NodeData())            # new root
        self.add(root)                              # add to tree description
        self.root = root.id                           # set as root

        # add branch to ingroup to unrooted tree
        self.unrooted.append(
            [root.id, ingroup_node, root_branch[2] - distance, root_branch[3]])
        # add branch to outgroup to unrooted tree
        self.unrooted.append([root.id, outgroup_node, distance, 0.0])
        _connect_subtree(root.id, ingroup_node)      # add ingroup
        _connect_subtree(root.id, outgroup_node)     # add outgroup

        # if theres still a lonely node in self.chain, then it's the old root,
        # and we delete it
        oldroot = [i for i in self.all_ids() if self.node(
            i).prev is None and i != self.root]
        if len(oldroot) > 1:
            raise ValueError('Isolated nodes in tree description: %s' %
                             ','.join(oldroot))
        elif len(oldroot) == 1:
            self.kill(oldroot[0])
        return self.root

    def to_string(self,
                  support_as_branchlengths=False,
                  branchlengths_only=False,
                  plain=True,
                  write_all_taxa=False,
                  branchlength_format="%1.5f",
                  support_format="%1.2f",
                  format="nexus"):
        """Return a paup compatible tree line.

        Arguments
        ---------
        support_as_branchlengths : bool
           If true, output bootstrap support value as branch lengths.
        branchlengths_only : bool
           Only output branchlengths, no support values
        plain : bool
           Output plain tree (no branch lengths/support values).
        write_all_taxa : bool
           If true, internal node names are output
        branchlength_format : string
           Format to use for branch lengths
        support_format : string
           Format to use for bootstrap support values
        format : string
           Either ``nexus`` on ``NHX``.

        Returns
        -------
        tree : string
            A PAUP compatible tree line.

        """
        # if there's a conflict in the arguments, we override plain=True
        if support_as_branchlengths or branchlengths_only:
            plain = False
        self.support_as_branchlengths = support_as_branchlengths
        self.branchlengths_only = branchlengths_only
        self.plain = plain
        self.write_all_taxa = write_all_taxa

        def make_info_string(data, terminal=False):
            """Creates nicely formatted support/branchlengths."""
            # CHECK FORMATTING

            if format == "nhx":

                info = "[&&NHX"
                try:
                    info += ":S=%s" % data.species
                except AttributeError:
                    pass

                if not terminal:
                    info += ":B=%s" % (support_format % data.support)
                return ':%s%s]' % (branchlength_format % data.branchlength,
                                   info)

            if self.plain:  # plain tree only. That's easy.
                return ''
            # support as branchlengths (eg. PAUP), ignore actual branchlengths
            elif self.support_as_branchlengths:
                if terminal:    # terminal branches have 100% support
                    return ':%1.2f' % (support_format % self.max_support)
                else:
                    return ':%1.2f' % (support_format % data.support)
            # write only branchlengths, ignore support
            elif self.branchlengths_only:
                return ':%s' % (branchlength_format % data.branchlength)
            # write suport and branchlengths (e.g. .con tree of mrbayes)
            else:
                if format in ("nh", "nexus"):
                    if terminal:
                        return ':%s' % (branchlength_format %
                                        data.branchlength)
                    else:
                        # we have blen and suppport
                        if data.branchlength is not None and \
                           data.support is not None:
                            sup = support_format % data.support
                            bl = branchlength_format % data.branchlength
                        # we have only blen
                        elif data.branchlength is not None:
                            sup = support_format % 0
                            bl = branchlength_format % data.branchlength
                        # we have only support
                        elif data.support is not None:
                            sup = support_format % data.support
                            bl = branchlength_format % 0
                        else:
                            sup = support_format % 0
                            bl = branchlength_format % 0
                        return "%s:%s" % (sup, bl)

        def newickize(node):
            """convert a node tree to a newick tree recursively.
            """
            if not self.node(node).succ:  # terminal
                return self.node(node).data.taxon + \
                    make_info_string(self.node(node).data, terminal=True)
            else:
                return '(%s)%s' % (','.join(map(newickize,
                                                self.node(node).succ)),
                                   make_info_string(self.node(node).data))

        def newickize_all_taxa(node):
            if not self.node(node).succ:    # terminal
                return self.node(node).data.taxon + \
                    make_info_string(self.node(node).data, terminal=True)
            else:
                if self.node(node).data.taxon is not None:
                    # changed output: first taxon name, then branch length
                    return '(%s)%s%s' % (
                        ','.join(map(newickize_all_taxa,
                                     self.node(node).succ)),
                        self.node(node).data.taxon,
                        make_info_string(self.node(node).data))
                else:
                    return '(%s)%s' % (
                        ','.join(map(newickize_all_taxa,
                                     self.node(node).succ)),
                        make_info_string(self.node(node).data))

        treeline = ''
        if format == "nexus":
            treeline = 'tree '
            if self.name:
                treeline += self.name
            else:
                treeline += 'a_tree'

            treeline += ' = '
            if self.weight != 1:
                treeline += '[&W%s] ' % str(round(float(self.weight), 3))
            if self.rooted:
                treeline += '[&R] '

        if self.write_all_taxa:
            treeline += '(%s);' % ','.join(map(newickize_all_taxa,
                                               self.node(self.root).succ))
        else:
            treeline += '(%s);' % ','.join(map(newickize,
                                               self.node(self.root).succ))
        return treeline

    def get_nodes(self, node_id):
        """Return a list of nodes downwards from a node (self, node_id).

        The list includes the given node_id.
        """
        if node_id is None:
            node_id = self.root
        if node_id not in self.chain:
            raise ValueError('Unknown node_id: %d.' % node_id)
        list = [node_id]
        for succ in self.chain[node_id].succ:
            list.extend(self.get_nodes(succ))
        return list

    def get_leaves(self, node_id):
        """Return a list of leaf nodes downward from a node (self, node_id).
        """
        return [x for x in self.get_nodes(node_id) if self.node(x).succ == []]

    def root_midpoint(self):
        """perform midpoint rooting of tree.

        The root is placed at equal distance to the two leaves
        furthest apart in the tree (centroid of the tree).
        """

        nnodes = max(self.chain.keys()) + 1

        # distance to root of parent including
        # the branch length of node to parent
        map_N2Root = [0] * (nnodes)
        # distance for node to leaves excluding branch length of this node
        map_N2Leaves = [0] * (nnodes)
        # distance for node to other part of trees excluding its children
        map_N2Other = [0] * (nnodes)

        def dist2leaves(node_id):
            node = self.node(node_id)
            if node.succ:
                map_N2Leaves[node_id] = max(
                    [map_N2Leaves[s] + self.node(s).data.branchlength for s in node.succ])

        def dist2root(node_id):
            node = self.node(node_id)
            map_N2Root[node_id] = node.data.branchlength
            if node.prev is not None:
                map_N2Root[node_id] += map_N2Root[node.prev]

        # traverse tree and record distance to root and maximum distance to
        # leaf at each node
        self.dfs(self.root,
                 pre_function=dist2root,
                 post_function=dist2leaves)

        # update distance to root for each node with length to leaves on the
        # other part of the tree.
        s = set(self.node(self.root).succ)
        for x in s:
            d = max(
                [map_N2Leaves[i] + self.node(i).data.branchlength for i in s.difference(set((x,)))])
            for xx in self.get_nodes(x):
                map_N2Root[xx] += d

        # update distance to any other
        def dist2other(node_id):
            node = self.node(node_id)
            if node.prev is not None:
                p = self.node(node.prev)
                root_dist = map_N2Root[node.prev]
                s = set(p.succ)
                d = max(
                    [map_N2Leaves[i] + self.node(i).data.branchlength for i in s.difference(set((node_id,)))])
                max_d = max(
                    (d, root_dist, map_N2Other[node.prev] + p.data.branchlength))
                map_N2Other[node_id] = max_d

        self.dfs(self.root,
                 pre_function=dist2other)

        # note: need not treat root. Set d2leaves > d2root
        # so that it is never chosen selected.
        map_N2Leaves[self.root] = 1
        map_N2Root[self.root] = 0

        # find node with minimum balance
        best_node = None
        best_balance = None
        for x in range(0, nnodes):
            if x == self.root:
                continue
            node = self.node(x)
            balance = math.fabs(map_N2Other[x] - map_N2Leaves[x])
            # need to add a small quantity to node.data.branchlength
            # to avoid rounding errors for large trees with uniform
            # branchlengths
            if balance >= 0 and balance <= node.data.branchlength + 0.001:
                if best_balance is None or best_balance > balance:
                    best_balance = balance
                    best_node = x

        # in trees with all 0 distances, best_node will be None
        # in this case return a balanced tree
        if best_node is None:
            return self.root_balanced()

        # compute distance of current node to new node
        d1 = map_N2Leaves[best_node]

        distance = float(map_N2Other[
                         best_node] - map_N2Leaves[best_node] + self.node(best_node).data.branchlength) / 2.0

        balance = map_N2Other[best_node] - map_N2Leaves[best_node]
        distance = float(
            balance + self.node(best_node).data.branchlength) / 2.0

        # remove current root
        self.unroot()

        # reroot tree
        result = self.root_at_node(best_node, distance=distance)

    def getNumLeaves(self):
        """return list with number of leaves beyond each node
        """
        nnodes = max(self.chain.keys()) + 1

        counts = [0] * nnodes

        def count(node_id):
            s = self.node(node_id).succ
            if s == []:
                counts[node_id] = 1
            else:
                for n in s:
                    counts[node_id] += counts[n]

        self.dfs(self.root, post_function=count)
        return counts

    def root_balanced(self):
        """perform balanced rooting of tree.

        The root is placed such that the number of leaves
        on either side of the tree is equal (or at most
        different by 1).
        """
        nnodes = max(self.chain.keys()) + 1
        nleaves = len(self.get_terminals())

        # number of leaves
        map_num_leaves = self.getNumLeaves()

        # simply find a node whose number of leaves is nleaves / 2
        threshold = int(math.floor(nleaves / 2.0))
        best_node = None
        for x in range(nnodes):
            if map_num_leaves[x] == threshold:
                best_node = x
                break

        if best_node is None:
            raise "no best node found for %i leaves" % threshold

        # remove current root
        self.unroot()

        # reroot tree
        self.root_at_node(best_node)

    def dfs(self, node_id,
            pre_function=Nop,
            descend_condition=Nop,
            post_function=Nop):
        """dfs tree tree traversal starting at node_id.

        Apply functions pre_function at first and
        post_function at last visit of a node.
        """
        pre_function(node_id)
        for n in self.node(node_id).succ:
            if descend_condition(node_id):
                self.dfs(n, pre_function, descend_condition, post_function)
        post_function(node_id)

    def writeToFile(self, outfile, with_branchlengths=True, format="nh"):
        """write a tree to a file."""

        if format == "nh":
            outfile.write(">%s\n%s\n" % (
                self.name,
                self.to_string(branchlengths_only=with_branchlengths,
                               write_all_taxa=True,
                               format="nh")))

        elif format == "nhx":
            outfile.write(">%s\n%s\n" % (self.name,
                                         self.to_string(format="nhx")))

    def truncate(self, node_id, taxon=None, keep_node=None):
        """truncate tree at node_id.

        This function will not change any branch lengths.
        If keep is given, single child nodes will be collapsed
        until keep_node is reached.
        """
        succ = tuple(self.node(node_id).succ)
        if len(succ) == 0:
            # if this is a leaf, do the same as tree.prune, but
            # nod not change branch length
            prev = self.unlink(node_id)
            self.kill(node_id)
            if not prev == self.root and len(self.node(prev).succ) == 1 \
                    and prev != keep_node:
                succ = self.node(prev).succ[0]
                self.collapse(prev)
        else:
            for s in succ:
                self.unlink(s)
                self.kill(s)

        if taxon:
            self.node(node_id).data.taxon = taxon

    def relabel(self, map_old2new, warn=False):
        """relabel taxa in tree using the provided mapping.
        """
        for node_id, node in list(self.chain.items()):
            if not node.data.taxon:
                continue

            if node.data.taxon in map_old2new:
                node.data.taxon = map_old2new[node.data.taxon]
            else:
                if warn:
                    raise KeyError("taxon %s not in map" % node.data.taxon)

    def rescaleBranchLengths(self, value):
        """rescale branch length so that they sum up to value."""

        t = 0.0
        for node_id, node in list(self.chain.items()):
            t += node.data.branchlength
        if t == 0:
            return
        for node_id, node in list(self.chain.items()):
            node.data.branchlength /= t
        return
