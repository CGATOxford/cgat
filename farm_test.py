import random, os, time

queue="test.q"

datafile = "tmp.data"
scriptfile = "tmp.scriptfile"
resultfile = "tmp.out"
logfile = "tmp.log"

num_jobs = 100
num_parts = 20
max_wait = 4

outfile = open( datafile, "w")
outfile.write("id\twait\n")
for x in range(0,num_parts):
    outfile.write( "%i\t%i\n" % (x, random.randint( 0, max_wait ) ) )
outfile.close()

outfile = open( scriptfile, "w")
outfile.write("""import sys,time
logfile=open("my.log", "w")
print"id\twaited"

for line in sys.stdin:
   if line.startswith("id"): continue
   id, wait = map(int, line[:-1].split("\t"))
   time.sleep(wait)
   print "%i\twaited %i" % (id, wait)
   logfile.write("logging for %i\\n"% id)

logfile.close()
""")

outfile.close()

cmd = "farm.py --cluster-priority=-10 --cluster-queue=%s --cluster-num-jobs=%i --split-at-column=1 --input-header --output-header --log=%s python %s --log=my.log < %s > %s" %\
               (queue, num_jobs, logfile, scriptfile, datafile, resultfile)
print cmd
os.system( cmd )

