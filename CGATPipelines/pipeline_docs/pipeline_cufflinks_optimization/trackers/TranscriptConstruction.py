from CGATReport.Tracker import *
import GTF
import IOTools


class TranscriptConstruction(TrackerSQL):

    '''
    tracker for testing transcript construction
    '''

    pattern = "(.*)_accepted_chr19_class"
    slices = ["protein_coding", "lincrna"]

    def __call__(self, track, slice=None):

        if slice is None:
            complete = self.getValue(
                """SELECT COUNT(*) FROM %(track)s_accepted_chr19_class WHERE class = 'complete'""")
            fragments = self.getValue(
                """SELECT COUNT(*) FROM %(track)s_accepted_chr19_class WHERE class like '%%fragment%%' """)
            return int(float(complete) / (complete + fragments) * 100)
        elif slice == "protein_coding":
            complete = self.getValue(
                """SELECT COUNT(*) FROM %(track)s_accepted_chr19_class WHERE class = 'complete' AND source = 'protein_coding' """)
            fragments = self.getValue(
                """SELECT COUNT(*) FROM %(track)s_accepted_chr19_class WHERE class like '%%fragment%%' AND source = 'protein_coding' """)
            return int(float(complete) / (complete + fragments) * 100)


class TranscriptConstruction2(TrackerSQL):

    '''
    tracker for testing transcript construction
    '''

    pattern = "(.*)_accepted_chr19_class"

    def __call__(self, track, slice=None):

        complete_gene_ids = self.getValues(
            """SELECT match_gene_id FROM %(track)s_accepted_chr19_class WHERE class = 'complete'""")
        fragment_gene_ids = self.getValues(
            """SELECT match_gene_id FROM %(track)s_accepted_chr19_class WHERE class like '%%fragment%%' """)

        frag = []
        for gene_id in fragment_gene_ids:
            if gene_id not in complete_gene_ids:
                frag.append(gene_id)
        complete = len(complete_gene_ids)
        fragments = len(frag)
        return int(float(complete) / (complete + fragments) * 100)


class LincRNADetection(TrackerSQL):

    '''
    tracker for testing transcript construction
    '''

    pattern = "(.*)_accepted_chr19_class"

    def getReferenceLincRNA(self, reference_gtf):

        lincs = []
        for entry in GTF.iterator(IOTools.openFile(reference_gtf)):
            if entry.source == "lincRNA":
                if entry.gene_id not in lincs:
                    lincs.append(entry.gene_id)
        return len(lincs)

    def __call__(self, track, slice=None):

        lincs = []
        lincs_detected = self.getValues(
            " SELECT match_gene_id FROM %(track)s_accepted_chr19_class WHERE (class = 'complete' OR class like '%%fragment%%') AND source = 'lincRNA' ")
        for linc in lincs_detected:
            if linc not in lincs:
                lincs.append(linc)

    #    return len(lincs)
        return int(float(len(lincs)) / float(self.getReferenceLincRNA("reference.chr19.gtf.gz")))
