import random
import sys
import gzip
import CGAT.Pipeline as P
import CGAT.Experiment as E

def write_random_records(fqa, fqb, outfa, outfb, N):
    """ get N random headers from a fastq file without reading the
    whole thing into memory"""

    records = sum(1 for _ in gzip.open(fqa)) / 4
    rand_records = sorted([random.randint(0, records - 1) for _ in xrange(N)])

    suba, subb = gzip.open(outfa, "w"), gzip.open(outfb, "w")
    fha, fhb = gzip.open(fqa),  gzip.open(fqb)
    rec_no = - 1
    
    for rr in rand_records:

        while rec_no < rr:
            rec_no += 1       
            for i in range(4): fha.readline()
            for i in range(4): fhb.readline()
        for i in range(4):
            suba.write(fha.readline())
            subb.write(fhb.readline())
        rec_no += 1 # (thanks @anderwo)

    print >>sys.stderr, "wrote to %s, %s" % (suba.name, subb.name)

if __name__ == "__main__":
    assert int(sys.argv[5]), "not a valid number to subsample"
    write_random_records(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
