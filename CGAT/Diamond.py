import random


class Diamond(object):

    def __init__(self):
        self.qid = None
        self.gi = None
        self.ref = None
        self.ngaps = 0
        self.length = 0
        self.evalue = 0
        self.nmismatches = 0
        self.identity = 0
        self.score = 0

    def read(self,
             qid,
             gi,
             ref,
             ngaps,
             length,
             evalue,
             nmismatches,
             identity,
             score):
        '''
        read input
        '''
        self.qid, \
            self.gi, \
            self.ref, \
            self.ngaps, \
            self.length, \
            self.evalue, \
            self.nmismatches, \
            self.identity, \
            self.score = \
            qid, \
            gi, \
            ref, \
            ngaps, \
            length, \
            evalue, \
            nmismatches, \
            identity, \
            score

        return self


def alignment_iterator(alignment_file):
    '''
    return an alignment
    '''
    # replaced .readlines()
    for line in alignment_file:
        data = line[:-1].split("\t")
        ref = data[1].split("|")
        if len(ref) == 1:
            # set ref to reference and leave gi as None
            gi, ref = None, ref[0]
        else:
            gi, ref = ref[1], ref[3]
        read, ngaps, length, evalue, nmismatches, identity, score = data[
            0], data[5], data[3], data[10], data[4], data[2], data[11]
        yield Diamond().read(read,
                             gi,
                             ref,
                             ngaps,
                             length,
                             evalue,
                             nmismatches,
                             identity,
                             score)


def query_iterator(alignment_iterator):
    '''
    iterate over alignments by query
    return a list of alignments
    '''
    reads_found = set()
    last = None
    alignments = []
    for alignment in alignment_iterator:
        this = alignment.qid
        if this != last:
            if last:
                yield alignments
            alignments = []
            reads_found.add(alignment.qid)
        alignments.append(alignment)
        last = this
    if last:
        yield alignments


def best_alignment_iterator(query_iterator):
    '''
    return stats for best alignment - best is determined
    by highest bit score. If there are multiple alignments
    with the same score return a random one
    '''
    for alignments in query_iterator:
        scores = []
        for alignment in alignments:
            scores.append(alignment.score)
        best = max(scores)
        best_alignments = [x for x in alignments if x.score == best]
        if len(best_alignments) > 1:
            best_alignments = random.sample(best_alignments, 1)
        best_alignment = best_alignments[0]
        yield best_alignment
