"""
=====================
Mappability pipeline
=====================

Pipeline to count mappable bases in a given genome

"""
import sys
import os

from ruffus import *
import CGAT.Experiment as E
import CGAT.Pipeline as P


###################################################
# Pipeline configuration
###################################################
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False},
    only_import=__name__ != "__main__")

PARAMS = P.PARAMS

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


@files(os.path.join(PARAMS["gem_dir"], PARAMS["genome"] + ".gem"),
       PARAMS["genome"] + ".mappability")
def calculateMappability(infile, outfile):
    '''Calculate mappability using GEM '''
    index = P.snip(infile, ".gem")
    to_cluster = True
    window_size = PARAMS("gem_window_size")
    threads = PARAMS("gem_threads")
    mismatches = PARAMS("gem_mismatches")
    max_indel_length = PARAMS("gem_max_indel_length")
    statement = '''gem-mappability
    -t %(threads)s -m %(mismatches)s
    --max-indel-length %(max_indel_length)s
    -l %(window_size)s -I %(index)s -o %(outfile)s ''' % locals()
    P.run()

###################################################################


@transform(calculateMappability, suffix(".mappability"), ".mappability.count")
def countMappableBases(infile, outfile):
    '''Count mappable bases in genome'''
    statement = '''cat %(infile)s | tr -cd ! | wc -c > %(outfile)s''' % locals()
    P.run()

###################################################################


@transform(countMappableBases, suffix(".count"), ".count.load")
def loadMappableBases(infile, outfile):
    '''load count of mappable bases in genome'''
    header = "total_mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=total_mappable_bases
                      --header-names=%(header)s
                   > %(outfile)s '''
    P.run()

###################################################################


@transform(calculateMappability, suffix(".mappability"), ".split.log")
def splitMappabiliyFileByContig(infile, outfile):
    '''Count mappable bases in genome'''
    track = P.snip(os.path.basename(infile), ".mappability")
    statement = '''mkdir contigs; 
                   csplit -k -f contigs/contig %(infile)s '/^~[a-zA-Z]/' {100000} > %(outfile)s;
                   rm contigs/contig00;''' % locals()
    P.run()

###################################################################


@follows(splitMappabiliyFileByContig)
@merge("contigs/contig*", PARAMS["genome"] + "_mappability_per_contig.tsv")
def countMappableBasesPerContig(infiles, outfile):
    '''Count mappable bases for each contig'''
    for infile in infiles:
        statement = '''grep '~' %(infile)s | sed s/~//g >> %(outfile)s; cat %(infile)s | tr -cd ! | wc -c >> %(outfile)s'''
        P.run()

    statement = '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s;'''
    P.run()

###################################################################


@transform(countMappableBasesPerContig, suffix(".tsv"), ".tsv.load")
def loadMappableBasesPerContig(infile, outfile):
    '''load count of mappable bases per contig '''
    header = "contig,mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=mappable_bases_per_contig
                      --header-names=%(header)s
                   > %(outfile)s '''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(calculateMappability, countMappableBases,
         loadMappableBases, splitMappabiliyFileByContig,
         countMappableBasesPerContig, loadMappableBasesPerContig)
def full():
    '''Count mappable bases in genome'''
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
