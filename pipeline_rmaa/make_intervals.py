#!/home/belgardt/bin/python2.6

# make a file giving interval categories with coordinates with the following hierarchy:
# coding sequence
# untranslated exons of protein-coding loci
# all other Ensembl annotations
# all introns
# expanded (intergenic regions in Ensembl that are connected to annotated regions by Cufflinks models)
# intergenic

# takes the following files:
# Ensembl GTF file
# Combined Cufflinks file (all.combined.gtf) | add all exons that are in the igi superstructure with an Ensembl annotation

# and returns a tab-separated file with all categories for all chromosomes used
# chr1  1   2186    intergenic
# chr1  2187    9058    expanded
# etc.

import interval
from sys import argv
from commands import getstatusoutput

ens_gtf_name = argv[1]
cuff_name = argv[2]
ens_gtf = open(ens_gtf_name,'r')
cuff = open(cuff_name,'r')

# get all chrom names
(status, output) = getstatusoutput("cat %s | awk '{ print $1 }' | sort | uniq" % ens_gtf_name )
if not status == 0:
    print "ERROR, CANNOT GET CHROM NAMES"
    sys.exit(2)
chroms = output.rstrip('\n').split('\n')

# initiate all categories as dictionaries of interval sets (chrom names are keys)
#CDS = pc_exons = cuff_exons = cuff_introns = other_ens = genic = intergenic = dict.fromkeys(chroms, interval.IntervalSet([]))
CDS = dict.fromkeys(chroms, interval()); pc_exons = dict.fromkeys(chroms, interval()); cuff_exons = dict.fromkeys(chroms, interval())
cuff_introns = dict.fromkeys(chroms, interval()); other_ens = dict.fromkeys(chroms, interval()); genic = dict.fromkeys(chroms, interval())
intergenic = dict.fromkeys(chroms, interval()); intronic = dict.fromkeys(chroms, interval())
# CDS: set of all protein-coding sequence
# pc_exons: set of all protein-coding exons (subtract CDS to get UTRs)
# other_ens: other Ensembl annotations
# genic: set of all genic regions - remainder after subtracting all exonic regions are intronic regions
# intronic: intronic regions
# cuff_introns: genic regions (intronic) which are expanded by Cufflinks models
# cuff_exons: genic regions (exonic) which are expanded by Cufflinks models
# intergenic: all intergenic regions
# UTRs: UTRs of pc genes (created later)

# iterate through the file grabbing CDS, pc sequence, other_ens, genic regions
gene_ids={}
for line in ens_gtf:
    la = line.rstrip('\n').split('\t')
    if la[1] == "protein_coding":
        if la[2] == "CDS":
            CDS[ la[0] ] = CDS.get(la[0], interval()) | interval[min( map(int,la[3:5]) ), max( map(int,la[3:5]) )]
        else:
            pc_exons[ la[0] ] = pc_exons.get(la[0], interval()) | interval[min( map(int,la[3:5]) ), max( map(int,la[3:5]) )]
    else:
        other_ens[ la[0] ] = other_ens.get(la[0], interval()) | interval[min( map(int,la[3:5]) ), max( map(int,la[3:5]) )]
    gene_id = la[8].split('";')[0].split('"')[1]
    gene_ids[ gene_id ] = gene_ids.get(gene_id, [la[0],set([])])
    gene_ids[ gene_id ][1].add(int(la[3])); gene_ids[ gene_id ][1].add(int(la[4]))

for gene_id, coords in gene_ids.iteritems():
    genic[coords[0]] = genic.get(coords[0], interval()) | interval[min(coords[1]),max(coords[1])]

# get all intronic
for chrom in chroms:
    # iterate through all intervals in genic[chrom], removing those pieces found in pc_exons and other_ens
    for interv in genic.get(chrom, interval()).components:
        if len( interv & pc_exons.get(chrom, interval()) ) == 0 and len( interv & other_ens.get(chrom, interval()) ) == 0:  # if no overlap at all (if we're good)
            intronic[chrom] = intronic.get(chrom, interval()) | genic[chrom]
        elif len( interv & pc_exons.get(chrom, interval()) ) == 0 and len( interv & other_ens.get(chrom, interval()) ) > 0: # if only overlaps other_ens
        elif len( interv & pc_exons.get(chrom, interval()) ) > 0 and len( interv & other_ens.get(chrom, interval()) ) == 0: # if only overlaps 
        else:   # if overlaps both
    intronic### NOTE STILL UPDATING FROM HERE DOWNWARD!!!!!!!!!
    intronic[chrom].difference_update( pc_exons[chrom] )
    intronic[chrom].difference_update( other_ens[chrom] )

UTRs = pc_exons
for chrom in chroms:
    UTRs[chrom].difference_update( CDS[chrom] )

old_idi_id=""; cuff_coords=interval(); firsttime=True; cuff_present=False
for line in cuff:
    la = line.rstrip('\n').split('\t')
    idi_id = la[8].split('";')[0].split('"')[1]
    if (not firsttime) and (not idi_id == old_idi_id):
        old_idi_id = idi_id
        if ens_here and cuff_present:
            cuff_exons[chrom].update( cuff_coords )
            cuff_introns_here = interval.IntervalSet([interval.Interval( cuff_exons[chrom].lower_bound(), cuff_exons[chrom].upper_bound() )])
            cuff_introns_here.difference_update( cuff_exons[chrom] )
            cuff_introns[chrom].update(cuff_introns_here)
        cuff_coords=interval.IntervalSet([])
        ens_here = False
        cuff_present = False
    firsttime = False
    chrom = la[0]
    if la[1] == "Cufflinks":
        cuff_coords.add( interval.Interval(min( map(int,la[3:5]) ), max( map(int,la[3:5]) )) )
        cuff_present=True
    else:
        ens_here = True
if ens_here and cuff_present:   # accounting for the final iteration
    cuff_exons[chrom].update( cuff_coords )
    cuff_introns_here = interval.IntervalSet([interval.Interval( cuff_exons[chrom].lower_bound(), cuff_exons[chrom].upper_bound() )])
    cuff_introns_here.difference_update( cuff_exons[chrom] )
    cuff_introns[chrom].update(cuff_introns_here)
for chrom in chroms:
    cuff_introns[chrom].difference_update(genic[chrom])
    cuff_exons[chrom].difference_update(genic[chrom])

# function to return a list of tuples (start, end, category name)
def return_intervals(interv_set, catname):
    coords = []
    for interv in interv_set:
        coords.append( (interv.lower_bound, interv.upper_bound, catname) )
    return coords

# CDS
# UTRs
# other_ens
# intronic
# cuff_exons
# cuff_introns
for chrom in chroms:
    sorted_coords = []
    sorted_coords.extend( return_intervals(CDS[chrom],"CDS") )
    sorted_coords.extend( return_intervals(UTRs[chrom],"UTR") )
    sorted_coords.extend( return_intervals(other_ens[chrom],"other_ens") )
    sorted_coords.extend( return_intervals(intronic[chrom],"intron") )
    sorted_coords.extend( return_intervals(cuff_exons[chrom],"cuff_exons") )
    sorted_coords.extend( return_intervals(cuff_introns[chrom],"cuff_introns") )
    sorted_coords.sort()
    for trio in sorted_coords:
        print "\t".join([chrom, str(trio[0]), str(trio[1]), trio[2]])

