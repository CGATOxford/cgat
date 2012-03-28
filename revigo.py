
import itertools, math, sys, optparse, collections

import GO, IOTools

import Experiment as E
import gographer.GOGraph
import numpy
import networkx.algorithms.traversal.depth_first_search as traversal

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

def outputList( s, go2info, p, lp ):

    for x in s: 
        print x, go2info.get(x, "na"), p.get(x, "na"), lp.get(x,"na")
        
def readSimrelMatrix( infile ):
    '''read simrel matrix from file.'''
    matrix, row_headers, col_headers = IOTools.readMatrix( infile )
    assert row_headers == col_headers
    return matrix, row_headers

def computeAncestors( graph, nodes ):
    '''build dictionary with node ancestors in graph.'''

    E.info("computing ancestral information")

    iteration, total = 0, len(nodes)

    ancestors = {}
    for node in nodes:
        a = getAncestors( graph, node )
        # only take nodes with annotations
        a = set([ x for x in a if x in nodes ])
        # add self
        # a.add( node )
        ancestors[node] = a 
        iteration += 1
        if iteration % 1000 == 0: 
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

def buildSimrelMatrix( graph, go2genes, go2info, ancestors, p ):
    '''build simrel matrix from an ontology DAG in *graph*
    and a mapping of go terms to genes.

    This method follows the procedure and syntax
    according to Schlicker et al. (2006) 
    BMC Bioinformatics 7:302
    '''

    ############################################################
    # find root
    root = getRoot( graph )

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

    with E.openOutputFile( "matrix" ) as outfile:
        IOTools.writeMatrix( outfile, matrix, test_nodes, test_nodes )

    return matrix, test_nodes

def main( argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"])

    parser.add_option("-o", "--filename-obo", dest="filename_obo", type="string",
                      help="filename with ontology information (in obo-xml format) [default=%default]." )

    parser.add_option("-g", "--filename-go", dest="filename_go", type="string",
                      help="filename with gene to go assignments [default=%default]." )

    parser.add_option("-m", "--filename-matrix", dest="filename_matrix", type="string",
                      help="filename with matrix [default=%default]." )

    parser.add_option( "--ontology", dest="ontology", type="choice", action="append",
                       choices=("biol_process","cell_location","mol_function", "mgi" ),
                       help="go ontologies to analyze. Ontologies are tested separately."
                       " [default=%default]." )

    parser.set_defaults( filename_obo = 'go_daily-termdb.obo-xml',
                         filename_go = "/ifs/data/annotations/mm9_ensembl64/go.tsv.gz",
                         #filename_pvalues = "/ifs/projects/proj008/report/genelists.go.dir/rela_upstream-changed.biol_process.results", 
                         filename_pvalues = "test.list", 
                         filename_bg = None, 
                         ontology = [],
                         # 0.9, 0.7, 0.5 and 0.4
                         # 0.4 very small clusters
                         min_simrel_threshold = 0.7,
                         min_frequency = 0.05,
                         min_pvalue = 0.05,
                         min_child_threshold = 0.75,
                         )

    options, args = E.Start( parser, argv, add_output_options = True )

    map_ontology2namespace = {
        'biol_process' :"biological_process",
        'cell_location': 'cellular_component',
        'mol_function': 'molecular_function' }

    for test_ontology in options.ontology:

        namespace = map_ontology2namespace[test_ontology]

        #########################################
        E.info("reading gene to go assignments from %s " % options.filename_go)
        all_gene2gos, all_go2infos = GO.ReadGene2GOFromFile( IOTools.openFile( options.filename_go ) )

        go2info = all_go2infos[test_ontology]
        gene2gos = all_gene2gos[test_ontology]

        #########################################
        with IOTools.openFile( options.filename_pvalues ) as infile:
            term2pvalue = {}
            # for line in infile:
            #     if line.startswith("code"): continue
            #     data = line[:-1].split()
            #     if len(data) < 14: continue
            #     pvalue = float( data[8] )
            #     if pvalue > options.min_pvalue: continue
            #     term2pvalue[data[11]] = pvalue
            for line in infile:
                if line.startswith("goid"): continue
                data = line[:-1].split()
                term2pvalue[data[0]] = float(data[1])
                
        if options.filename_bg:
            background = IOTools.readList( IOTools.openFile( filename_bg )) 
        else:
            background = list(gene2gos.keys())

        E.info( "reading GO graph from %s" % options.filename_obo)
        graph = gographer.GOGraph( namespace = namespace, GOOboXmlFileName=options.filename_obo )

        # compute ancestry information for all nodes
        ancestors = computeAncestors( graph, graph.nodes() )
        
        # invert mapping, at the same time insert missing ancestral counts
        go2genes = GO.buildGO2Genes( gene2gos, ancestors = ancestors )
        counts, p = computeTermFrequencies( graph, go2genes )

        if options.filename_matrix:
            with IOTools.openFile( options.filename_matrix, "r"):
                matrix, terms = readSimrelMatrix( infile )
        else:
            matrix, terms = buildSimrelMatrix( graph, go2genes, go2info, ancestors, p )

        E.info("working on simrel matrix with %i terms" % len(terms))
        term2index = dict( [ (y,x) for x,y in enumerate( terms ) ] )

        # shrink matrix to only those terms that are in result set
        terms = [ x for x in terms if x in term2pvalue ]

        E.info("after removing terms with no p-value: %i terms" % len(terms))
        if len(terms) == 0:
            raise ValueError( "no terms given after p-value filtering" )
        matrix = matrix[ numpy.array( [ term2index[x] for x in terms ]) ]
        matrix = matrix[ :, numpy.array( [ term2index[x] for x in terms ]) ]
        term2index = dict( [ (y,x) for x,y in enumerate( terms ) ] )

        # cluster according to revigo
        clusters = dict( [ (x, [x]) for x in terms ] )
        nrows = len(terms)
        
        counter = E.Counter()

        print matrix[term2index['GO:0012507']][term2index['GO:0030134']]

        while 1:
            # find most similar pair of GO terms ti and tj with similarity simrel
            index = numpy.argmax( matrix )
            i = index // nrows
            j = index % nrows
            
            simrel = matrix[i][j]
            
            t1 = terms[i]
            t2 = terms[j]
            print "---------------------"
            print index, i, j, simrel
            assert simrel == matrix.max()

            print t1, go2info.get(t1, "na"), term2pvalue[t1], p[t1] * 100.0
            print t2, go2info.get(t2, "na"), term2pvalue[t2], p[t2] * 100.0

            p1 = term2pvalue[t1]
            p2 = term2pvalue[t2]

            # stop procedure if similarity is below threshold
            if simrel < options.min_simrel_threshold: 
                break

            counter.tested += 1
        
            if p[t1] > options.min_frequency:
                counter.min_frequency += 1
                # remove term with frequency > 5%
                rep, mem = t2, t1
            elif p[t2] > options.min_frequency:
                counter.min_frequency += 1
                # remove term with frequency > 5%
                rep, mem = t1, t2
            elif term2pvalue[t1] < term2pvalue[t2]:
                counter.min_pvalue += 1
                # remove term with higher pvalue
                rep, mem = t1, t2
            elif term2pvalue[t2] < term2pvalue[t1]:
                counter.min_pvalue += 1
                # remove term with higher pvalue
                rep, mem = t2, t1
            elif t1 in ancestors[t2]:
                # reject child term
                # exception: if parent is > 75% child, keep child
                ovl = float(len(counts[t1].intersection( counts[t2] )))
                if (ovl / len(counts[t1])) > options.min_child_threshold:
                    counter.ancestor_child += 1
                    mem, rep = t2, t1
                else:
                    counter.ancestor_parent += 1
                    rep, mem = t1, t2
            elif t2 in ancestors[t1]:
                # reject child term
                # exception: if parent is > 75% child, keep child
                ovl = float(len(counts[t1].intersection( counts[t2] )))
                if (ovl / len(counts[t2])) > options.min_child_threshold:
                    counter.ancestor_child += 1
                    mem, rep = t1, t2
                else:
                    counter.ancestor_parent += 1
                    rep, mem = t2, t1
            else:
                counter.random += 1
                # keep random term
                rep, mem = min(t1,t2), max(t1,t2)

            # merge clusters
            clusters[rep].extend( clusters[mem] )
            del clusters[mem]
            # set values of mem to 0
            if mem == t1: 
                matrix[i,:] = 0
                matrix[:,i] = 0
            else:
                matrix[j,:] = 0
                matrix[:,j] = 0

        E.info( "%s" % str(counter))
        E.info( "after clustering: %i clusters" % len(clusters) )

        for rep, members in clusters.iteritems():
            print rep, p[rep], go2info.get(rep, "na"), term2pvalue[rep]
            print "-->", members

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
