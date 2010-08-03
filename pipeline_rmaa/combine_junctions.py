#!/home/belgardt/bin/python2.6

from sys import stdin

for line in stdin:
    la=line.rstrip('\n').split('\t')
    my_str = la[0] + '\t' + \
        str( int(la[1]) + int( la[10].split(',')[0] ) - 1 ) + '\t' + \
        str( int( la[2] ) - int( la[10].split(',')[1] ) + 1 ) + '\t' + la[5]
    print my_str
