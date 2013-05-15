################################################################################
#   Gene prediction pipeline 
#
#   $Id: count_terminal_exons.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""quick script to count terminal exons.

input:

file1: list if with number of exons per token
file2: list of token exon id.

"""

if __name__ == "__main__":
    
    infile_num_exons = "num_exons_per_reference"
    infile_missed_exons = "missed_exon_ids"
    
    num_exons = {}

    for line in open(infile_num_exons, "r").readlines():
        if line[0] == "#": continue
        x,y = line.split("\t")[:2]
        num_exons[x] = int(y)

    missed_exons = {}
    
    for line in open(infile_missed_exons, "r").readlines():
        if line[0] == "#": continue
        x,y = line.split("\t")[:2]
        if x not in missed_exons: missed_exons[x] = []
        missed_exons[x].append(int(y))

    total_nterminal = 0
    total_cterminal = 0
    total_internal = 0
    total_missed = 0
    total = 0
    for x,yy in missed_exons.items():
        yy.sort()

        nnterminal = 0
        ncterminal = 0
        ninternal = 0
        
        total += len(yy)
        
        print "# %s\tnterminal=" % x ,
        index = 1
        for a in yy:
            if a == index:
                print "%i" % a,
                nnterminal += 1
                index += 1
            else:
                break
        index_n = index
        
        print "\tcterminal=",
        index = num_exons[x]
        yy.reverse()
        for a in yy:
            if a == index:
                print "%i" % a,
                ncterminal += 1
                index -= 1
            else:
                break
        index_c = index
        
        print "\tinternal=",
        yy.reverse()
        for a in yy[nnterminal:len(yy) - ncterminal]:
            print "%i" % a,
            ninternal += 1
        print
            
        status = "all"
        
        if nnterminal == num_exons[x]:
            status = "missed"
            total_missed += nnterminal
        else:
            if ninternal == 0 and nnterminal == 0 and ncterminal > 0 :
                status = "only-c"
            elif ninternal == 0 and nnterminal > 0 and ncterminal == 0:
                status = "only-n"
            elif nnterminal == 0 and ncterminal == 0:
                status = "only-i"
            
            total_cterminal += ncterminal
            total_nterminal += nnterminal            
            total_internal  += ninternal
            
        print "%s\t%s\t%i\t%i\t%i\t%i" % (x, status, num_exons[x], nnterminal, ncterminal, ninternal)

    print "total\tmissed\tnterminal\tcterminal\tinternal"
    print "%i\t%i\t%i\t%i\t%i" % (total, total_missed, total_nterminal, total_cterminal, total_internal)
    print "%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f" % \
          tuple(map( lambda x: 100 * float(x) / total,
                     (total, total_missed, total_nterminal, total_cterminal, total_internal) ))



