"""Class to compute connected components

Connected components list all nodes in a graph that
are connected either directly or via intermediate egdes.

The module computes connected components using the
algorithm listed by Sedgwick: Algorithms in C, p507.
Components are computed in an on-line fashion, so that
the whole graph is not kept in memory.

USAGE:

# create object for string identifiers
>>> x = SComponents()

# add edges
>>> x.add( "1", "2" )
>>> x.add( "1", "3" )
>>> x.add( "4", "5" )

# retrieve components
>>> print x.getComponents()

This is a cython extension class."""

cdef extern from "connected_components.h":

    ctypedef struct cSComponents "CharComponents":
        int add(char *, char *)
        int get(char *)
        int getComponent(int)
        int getIndex(char *)
        char * getToken(int)
        int getNumNodes()
        void reset()

    cSComponents *new_SComponents "new CharComponents" ()
    void del_SComponents "delete" (cSComponents * c)

    ctypedef struct cIComponents "IntComponents":
        int add(int, int)
        int get(int)
        int getComponent(int)
        int getIndex(int)
        int getToken(int)
        int getNumNodes()
        void reset()

    cIComponents *new_IComponents "new IntComponents" ()
    void del_IComponents "delete" (cIComponents * c)

cdef class SComponents:
    """Components for a graph with string indices."""
    cdef cSComponents *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self ):
        self.thisptr = new_SComponents()
    def __dealloc__(self):
        del_SComponents(self.thisptr)
    def add(self, a, b):
        """add an edge between nodes a and b
        return True, if the link joins two previously disconnected componenents.
        """
        return self.thisptr.add(a, b) > 0
    def getNumNodes(self):
        """return the number of nodes in the graph."""
        return self.thisptr.getNumNodes()
    def getComponent(self, id):
        """return component an id belongs to.
        If the token is unknown, return 0.
        """
        return self.thisptr.getComponent(id)
    def getToken(self, id):
        """return token for an id.
        """
        return self.thisptr.getToken(id)
    def getIndex(self, token):
        """return id for a token.
        If the token is unknown, return 0.
        """
        return self.thisptr.getIndex(token)
    def get(self, v):
        """return component for a token.
        If the token is unknown, return 0.
        """
        return self.thisptr.get(v)

    def getComponents(self):
        """return all connected components as a list of lists."""
        components = {}
        
        for x in range(1, self.thisptr.getNumNodes() + 1):
            c = self.thisptr.getComponent(x)
            if c not in components:
                components[c] = []
            components[c].append(self.thisptr.getToken(x))
        
        return components.values()

    def reset(self):
        """clear graph.
        """
        self.thisptr.reset()

cdef class IComponents:
    """Components for a graph with numeric indices."""
    cdef cIComponents *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self ):
        self.thisptr = new_IComponents()
    def __dealloc__(self):
        del_IComponents(self.thisptr)
    def add(self, a,b):
        """add an edge between nodes a and b
        return True, if the link joins two previously disconnected componenents.
        """
        return self.thisptr.add( a, b ) > 0
    def getNumNodes(self):
        """return the number of nodes in the graph."""
        return self.thisptr.getNumNodes()
    def getComponent(self, id):
        """return component an id belongs to.
        If the token is unknown, return 0.
        """
        return self.thisptr.getComponent(id)
    def getToken(self, id):
        """return token for an id.
        """
        return self.thisptr.getToken(id)
    def getIndex(self, token):
        """return id for a token.
        If the token is unknown, return 0.
        """
        return self.thisptr.getIndex(token)
    def get(self, v):
        """return component for a token.
        If the token is unknown, return 0.
        """
        return self.thisptr.get(v)

    def getComponents(self):
        """return all connected components as a list of lists."""
        components = {}
        
        for x in range(1, self.thisptr.getNumNodes() + 1):
            c = self.thisptr.getComponent(x)
            if c not in components: components[c] = []
            components[ c  ].append( self.thisptr.getToken( x ) )
        
        return components.values()

    def reset(self):
        """clear graph.
        """
        self.thisptr.reset()
