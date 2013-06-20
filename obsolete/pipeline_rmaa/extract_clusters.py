#!/home/belgardt/bin/python2.6
# extracts clusters from fuzzyK results

from sys import argv

result_file = open( argv[1], 'r' )  # output of FuzzyK
output_prefix = argv[2]         # common prefix for the output files
min_identity = float( argv[3] ) # min fraction identity required to cluster

num_centroids = len(result_file.readline().rstrip('\n').split('\t')) - 3  # get rid of header
outfiles = []
for centroid in range(num_centroids):
    outfiles.append( open( argv[2] + '.' + str(centroid) , 'w' ) )
for line in result_file:
    la = line.rstrip('\n').split('\t')
    gene = la[0]; values = la[3::]
    for centroid, value in enumerate(values):
        if float(value) >= min_identity:
            outfiles[ centroid ].write( gene + '\n' )

