################################################################################
#   Gene prediction pipeline 
#
#   $Id: Tree_test.py 2784 2009-09-10 11:41:14Z andreas $
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
"""unit testing module for the Tree.py class."""

from Tree import Tree
import unittest

class MidPointRootingCheck(unittest.TestCase):

    trees = [
        "(A:1,((B:1,C:1):1,D:1):1,(E:1,F:1):1);",
        "((A:1,(B:1,C:1):1):0.5,((D:1,E:1):1,F:1):0.5);",
        "(A:1,B:1);",
        "((A:1,B:1):1,C:5);",
        "((A:1,B:1):10,(D:1,C:5):1);",
        "((A:0,B:0):1,(D:0,C:5):0);",
        "((A:5,B:1):1,D:1);",
        "((A:0,(B:0,(((C:0,D:0):0,E:2.11270):0,((F:0,G:0):0,(H:0,I:0):0):0):0):0):0,((J:0,K:0.12496):0,L:0):0,M:0);",
        "(ID000001:1.19640,(((ID000004:0.41850,ID000006:0.06490):0.12010,ID000005:0.30820):0.31570,ID000003:0.38540):0.00000,ID000002:1.27200);",
        "(A:0.19174,(B:0.58034,((C:0.98961,D:0.52099):0.14598,(E:0.00000,F:0.00000):0.49107):0.67347):0.01248,G:0.34146);",
        "(((((A:0.10670,B:0.35050):0.03480,C:0.13720):0.23850,D:0.31120):0.12570,E:0.38110):0.04180,F:0.79130);",
        "(ID000001:0.83310,(((ID000005:0.10670,ID000004:0.35050):0.03480,ID000006:0.13720):0.23850,ID000003:0.31120):0.12570,ID000002:0.38110);",
    ]
    
    def testLength(self):
        """midpoint rooting."""
        for tree in self.trees:
            
            t = Tree( tree )
            
            t.root_midpoint()
            # test 1: only two children for root
            s = t.node(t.root).succ
            self.assertEqual( len(s), 2 )

            # calculate tree length on either side of tree        
            d2leaves = [ 0 ] * (max( t.chain.keys() ) + 1)
            def dist2leaves( node_id ):
                node = t.node(node_id)
                if node.succ:
                    d2leaves[node_id] = max( [d2leaves[s] + t.node(s).data.branchlength for s in node.succ ] )

            t.dfs( t.root, post_function = dist2leaves )

            d1 = d2leaves[s[0]] + t.node(s[0]).data.branchlength
            d2 = d2leaves[s[1]] + t.node(s[1]).data.branchlength
            # test 2: distance to children equal on both sides
            self.assertAlmostEqual( d1, d2, 5,
                                    "assertion error: %s != %s for tree %s -> %s" % \
                                    (str(d1), str(d2),
                                     tree, t.to_string( branchlengths_only=True ) ))

    def testNegativeBranchLengths(self):
        for tree in self.trees:
            
            t = Tree( tree )

            t.root_midpoint()

            for n, node in t.chain.items():
                self.failIf( node.data.branchlength < 0,
                             "assertion error: negative branchlength for tree %s -> %s" %\
                             (tree, t.to_string( branchlengths_only = True )) )
        
    def testTruncate(self):
        
        t = Tree( "(A:1,((B:1,C:1):1,D:1):1,(E:1,F:1):1);" )
        t.truncate( 7, "E+F" )
        result = t.to_string(branchlengths_only = True,
                             branchlength_format = "%i",
                             format="nh")
        self.assertEqual( result, "(A:1,((B:1,C:1):1,D:1):1,E+F:1);" )

        t = Tree( "(A:1,((B:1,C:1):1,D:1):1,(E:1,F:1):1);" )
        t.truncate( 8 )
        result = t.to_string(branchlengths_only = True,
                             branchlength_format = "%i",
                             format="nh")
        self.assertEqual( result, "(A:1,((B:1,C:1):1,D:1):1,F:1);" )
        

if __name__ == "__main__":
        unittest.main()
