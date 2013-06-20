#!/home/belgardt/bin/python2.6

# Make a tree from a matrix of RPKM values weighing genes by transcription level

from sys import argv
from math import log
from commands import getoutput
from string import letters
from random import choice

infile = open(argv[1], 'r')
outfile = open(argv[2] + '.distance', 'w')

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
    for i in range(num_samples):
        for j in range(num_samples):
            the_matrix[i][j] += abs( la[i] - la[j] )

# write distance matrix
for i in range(num_samples):
    outfile.write( "%-10s" % header[i] )
    for j in range(num_samples):
        outfile.write( ' ' + str( the_matrix[i][j] ) )
    outfile.write( '\n' )
infile.close(); outfile.close()

# create tmp directory & work there - different syntax though if absolute vs relative path
# make commands file for fitch & run
commands_file = open( argv[2] + '.commands', 'w')
TMP_DIR = "".join([choice(letters) for x in xrange(10)]); getoutput('mkdir %s' % TMP_DIR)
if ABS_PATH:
    commands_file.write( '%s\nG\nJ\n23\n5000\n2\nY\n' % (argv[2] + '.distance') )
    commands_file.close()
    getoutput('cd %s; fitch < %s; rm outfile; mv outtree %s; cd ..' % ( TMP_DIR, argv[2] + '.commands', argv[2] ) )
else:
    commands_file.write( '../%s\nG\nJ\n23\n5000\n2\nY\n' % (argv[2] + '.distance') )
    commands_file.close()
    getoutput('cd %s; fitch < ../%s; rm outfile; mv outtree ../%s; cd ..' % ( TMP_DIR, argv[2] + '.commands', argv[2] ) )
getoutput('rmdir %s' % TMP_DIR )

