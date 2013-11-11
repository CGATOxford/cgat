"""
add def main() to scripts()
"""


import sys, re

scripts = sys.argv[1:]

TEXT1 = '''
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv
'''

TEXT2='''
if __name__ == "__main__":
    sys.exit( main( sys.argv) )
'''

for script in scripts:

    lines = open( script ).readlines()

    # check if def main( already exists
    if [ x for x in lines if x.startswith("def main(" ) ]:
        continue

    # find "if __name__"
    l = [ (x,y) for x,y in enumerate(lines) if y.startswith("if __name__" ) ]

    if len(l) == 0:
        print script, "contains no __name__"
        continue

    if len(l) > 1:
        print script, "contains multiple __name__"
        continue

    main_line = l[0][0]

    del lines[main_line]
    lines[main_line:main_line] = TEXT1
    lines.extend( TEXT2 )

    outf = open( script, "w" )
    outf.write( "".join(lines) + "\n" )
        
    print "updated: %s" % script

