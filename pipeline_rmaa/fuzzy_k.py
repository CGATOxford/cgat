'''setup an aerie run.'''

import Experiment as E

from sys import argv
from commands import getoutput
from os import path
from string import letters
from random import choice
from math import log
import sys, tempfile, subprocess, shutil, os, re

rpkm_filename, cdt_filename = argv[1], argv[2]
centroid_filename, membership_filename, background_filename = argv[3], argv[4], argv[5]

rpkm_file = open( rpkm_filename, 'r' )
cdt_file = open( cdt_filename, 'w')
centroid_file = argv[3]
membership_file = argv[4]
# write names of genes that meet the min length and min RPKM requirements, but not the min diff requirement
background_file = open( background_filename, 'w')   
min_gene_len = int(argv[6])
min_rpkm = float(argv[7])
min_diff = log(float(argv[8]),2)    # minimum fold difference between the highest and lowest to be considered (filter)

# make cdt file
labels = rpkm_file.readline().rstrip('\n').split('\t')[2::]
my_str = "UID\tNAME\tGWEIGHT"
my_str2 = "EWEIGHT\t\t"
for label in labels:
    my_str += "\t" + label
    my_str2 += "\t1"
cdt_file.write(my_str + '\n' + my_str2)

counts = E.Counter()

for line in rpkm_file:
    data = line[:-1].split("\t")
    counts.input += 1
    if int(data[1]) < min_gene_len:    # exclude genes that are too short
        counts.skipped_length += 1
        continue
    name = data[0]
    la = map(float, data[2:])

    if max(la) < min_rpkm: 
        counts.skipped_rpkm += 1
        continue

    background_file.write(name + "\n")

    la = map(lambda x: x + 0.01, la)    # to handle any zero values, add 0.01 to every RPKM
    avg_rpkm = float(sum(la) ) / len(la)
    ratios = [ log( x/avg_rpkm, 2) for x in la]

    if max(ratios) < min_diff:
        counts.skipped_diff += 1
        continue

    cdt_file.write("\n" + name + "\t" + name + "\t1\t" + "\t".join(map(str,ratios)) )

    counts.output += 1

cdt_file.close()

E.info( "counts: %s" % str(counts ) )
if counts.output == 0:
    E.warn("no genes output - fuzzy_k clustering aborted" )
    sys.exit(0)

# run fuzzy k on everything


nclusters = len(labels)+10
old_place = path.abspath('')

## build instruction file
my_str = '''load %s
fuzzy %s
alls
exit
''' % ( argv[2].rsplit('/', 1)[1], nclusters)

instruction_filename = os.path.abspath( 'fuzzy_k/fuzzykinstructions-%s-%s' % (str(min_rpkm), argv[8]) )

inst_file = open( instruction_filename,'w')
inst_file.write(my_str)
inst_file.close()

# run aerie in a temporary directory
tmpdir = tempfile.mkdtemp()

# move cdt_file locally (is this really necessary)
shutil.copy( argv[2], os.path.join( tmpdir, os.path.basename( argv[2] ) ) )

# Grant used to copy the excecutable into the local directory, too



try:
    retcode = subprocess.call( "arie < %(instruction_filename)s" % locals(), 
                               shell=True,
                               cwd = tmpdir )
    if retcode < 0:
        print >>sys.stderr, "Child was terminated by signal", -retcode
    else:
        print >>sys.stderr, "Child returned", retcode
except OSError, e:
    print >>sys.stderr, "Execution failed:", e
                           
shutil.move( os.path.join( tmpdir, "alls.fct"), centroid_filename )
shutil.move( os.path.join( tmpdir, "alls.mb"), membership_filename )
shutil.unlink( tmpdir )

