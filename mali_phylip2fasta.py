import sys, string, re

if __name__ == "__main__":

    line = sys.stdin.readline()
    num_lines, width = map(int, re.search("(\d+)\s+(\d+)", line[:-1]).groups())

    alignment = []

    for x in range(num_lines):
        line = sys.stdin.readline()
        id, ali = re.search( "^(\S+)\s+(.+)", line[:-1]).groups()
        alignment.append( (id, [ali]) )

    while 1:
        line = sys.stdin.readline()
        if not line: break
        for x in range(num_lines):
            line = sys.stdin.readline()            
            alignment[x][1].append( line[:-1] )

    for x in range(num_lines):
        print ">%s\n%s" % (alignment[x][0], re.sub("\s", "", string.join(alignment[x][1], "")))
