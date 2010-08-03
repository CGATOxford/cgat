#!/home/belgardt/bin/python2.6

# Calculates the overall total expression level of genes in each category
# Calculates the largest differences (enrichments and depletions) of categories in terms of overall expression (just compared to overall per gene level expression - NOT weighing genes equally)

from sys import argv

# 1: unique.rpkm.all
rpkm_file = open(argv[1],'r')
# 2: annotations file
annotations_file = open(argv[2],'r')
# 3: first output file - sorted by expression
expression_file = open(argv[3],'w')
# 4: second output file - sorted by enrichment (depletion at bottom)
diff_file = open(argv[4],'w')

# get average RPKMs for all genes
line=rpkm_file.readline(); la=line.rstrip('\n').split('\t'); num_samples=len(la[2::])
genes = {}
for line in rpkm_file:
    la = line.rstrip('\n').split('\t')
    genes[ la[0] ] = sum(map(float, la[2::]))/float(num_samples)

# get lists of average RPKMs for all annotations
annotations_file.readline()
annotations_rpkms = {}
for line in annotations_file:
    la = line.rstrip('\n').split('\t')
    if la[1] == "" or len(la)<2:
        continue
    annotations_rpkms[la[1]] = annotations_rpkms.get(la[1], [])
    try:
        annotations_rpkms[la[1]].append( genes[la[0]] )
    except:
        pass

# return annotations in descending order of total RPKM
annotation_tuples = []
for annotation, RPKMs in annotations_rpkms.iteritems():
    annotation_tuples.append( (sum(RPKMs), annotation) )
annotation_tuples.sort()
annotation_tuples.reverse()
for annotation in annotation_tuples:
    expression_file.write( str(annotation[0]) + '\t' + annotation[1] + '\n' )

# return annotations in descending order of enrichment
enrich_tuples = []
for annotation in annotation_tuples:
    if len(annotations_rpkms[annotation[1]]) == 0:
        continue
    enrich_tuples.append( (annotation[0]/len(annotations_rpkms[annotation[1]]), annotation[0], len(annotations_rpkms[annotation[1]]), annotation[1] ) )  # avg RPKM per gene, total RPKM, total # genes, term name
enrich_tuples.sort()
enrich_tuples.reverse()

for enriched in enrich_tuples:
    diff_file.write( "\t".join(map(str, enriched)) + '\n' )
