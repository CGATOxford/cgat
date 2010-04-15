import string, re, sys

if __name__ == "__main__":
    param_filename_patterns = sys.argv[1]

    infile = open(param_filename_patterns, "r")

    patterns = []
    for line in infile:
        (old, new) = line[:-1].split("\t")
        patterns.append( (re.compile( "%s" % old), new) )

    infile.close()

    for line in sys.stdin:

        for pattern, new in patterns:
            line = pattern.sub( new, line )

        print line[:-1]
