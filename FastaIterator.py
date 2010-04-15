"""iterate over fasta files.

The difference to the biopython iterator is that this one skips over comment lines
starting with "#".
"""

class FastaRecord:

    def __init__(self, title, sequence ):

        self.title = title
        self.sequence = sequence

def iterate( infile ):

    h = infile.readline()[:-1]

    if not h: raise StopIteration
    
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h: raise StopIteration
        continue

    h = h[1:]
    
    seq = []
    
    for line in infile:

        if line[0] == "#": continue
        
        line = line[:-1] # remove newline
        if not line: continue

        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            
            h = line[1:]
            seq = []
            continue
        
        seq.append(line)
        
    yield FastaRecord(h,''.join(seq))

class FastaIterator:

    def __init__(self, f, *args, **kwargs):
        self.mIterator = iterate(f)
        
    def __iter__(self):
        return self

    def next(self):
        return self.mIterator.next()

##------------------------------------------------------------
def iterate_together( *args ):
    """iterate over one or more fasta files and output tuples of sequences."""
    iterators = [ FastaIterator( open(x, "r") ) for x in args ]

    while 1:
        yield [ x.next() for x in iterators ]
