import ncl

class IndexedGenome:
    '''Genome with indexed intervals.
    '''

    index_factory = ncl.NCL
    def __init__( self ): 
        self.mIndex = {}

    def add( self, contig, start, end, value ):
        
        if contig not in self.mIndex: 
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add( start,end,value )

    def __getitem__(self, args):
        '''return intervals overlapping with key.'''
        if args[0] not in self.mIndex: 
            raise KeyError("contig %s not in index" % args[0] )
        
        return self.mIndex[args[0]].find( args[1],args[2] )

    def get(self, contig, start, end):
        '''return intervals overlapping with key.'''
        if contig not in self.mIndex: 
            raise KeyError("contig %s not in index" % contig )
        
        return self.mIndex[contig].find( start, end )

    def __len__(self):
        '''return number of contigs.'''
        return len(self.mIndex)

class Simple( IndexedGenome ):
    '''index intervals without storing a value.'''
    index_factory = ncl.NCLSimple
    
    def __init__(self, *args, **kwargs ):
        IndexedGenome.__init__(self, *args, **kwargs )

    def add( self, contig, start, end ):
        
        if contig not in self.mIndex: 
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add( start,end )
