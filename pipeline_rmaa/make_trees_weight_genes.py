#!/home/belgardt/bin/python2.6

# Make a tree from a matrix of RPKM values giving all genes equal weight (requiring at least one lane to be above a certain minimum RPKM cutoff)

from sys import argv
from math import log
from commands import getoutput
from string import letters
from random import choice

infile = open(argv[1], 'r')
min_rpkm = float(argv[2])
outfile = open(argv[3] + '.distance', 'w')

# try to detect if relative or absolute path
if argv[1][0] == '/' or argv[1][0] == '~':
    ABS_PATH = True
else:
    ABS_PATH = False

# write header
header = infile.readline().rstrip('\n').split('\t')[2::]
num_samples=len(header)
outfile.write("   %s\n" % num_samples)

# initialize output matrix
the_matrix=[]
for i in range(num_samples):
    the_matrix.append([0.0]*num_samples)

# build output matrix
for line in infile:
    la = map(float, line.rstrip('\n').split('\t')[2::])
    if max(la) < min_rpkm:
        continue
    la = map(lambda x: x + 0.01, la)    # to handle any zero values, add 0.01 to every RPKM
    avg_rpkm = float(sum(la))/len(la)
    ratios = map(lambda x: log(x/avg_rpkm, 2), la)
    for i in range(num_samples):
        for j in range(num_samples):
            the_matrix[i][j] += abs( ratios[i] - ratios[j] )

# write distance matrix
for i in range(num_samples):
    outfile.write( "%-10s" % header[i] )
    for j in range(num_samples):
        outfile.write( ' ' + str( the_matrix[i][j] ) )
    outfile.write( '\n' )
infile.close(); outfile.close()

# create tmp directory & work there - different syntax though if absolute vs relative path
# make commands file for fitch & run
commands_file = open( argv[3] + '.commands', 'w')
TMP_DIR = "".join([choice(letters) for x in xrange(10)]); getoutput('mkdir %s' % TMP_DIR)
if ABS_PATH:
    commands_file.write( '%s\nG\nJ\n23\n5000\nP\n0\n2\nY\n' % (argv[3] + '.distance') )
    commands_file.close()
    getoutput('cd %s; fitch < %s; rm outfile; mv outtree %s; cd ..' % ( TMP_DIR, argv[3] + '.commands', argv[3] ) )
else:
    commands_file.write( '../%s\nG\nJ\n23\n5000\nP\n0\n2\nY\n' % (argv[3] + '.distance') )
    commands_file.close()
    getoutput('cd %s; fitch < ../%s; rm outfile; mv outtree ../%s; cd ..' % ( TMP_DIR, argv[3] + '.commands', argv[3] ) )
getoutput('rmdir %s' % TMP_DIR )

