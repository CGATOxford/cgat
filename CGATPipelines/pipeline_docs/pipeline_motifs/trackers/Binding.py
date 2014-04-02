from IntervalReport import *


class BindingPatterns(IntervalTracker):

    '''output summary counts of binding patterns.

    The empty pattern is excluded.
    '''
    pattern = "(.*)_binding"

    def __call__(self, track):
        return self.getAll( """SELECT pattern, COUNT(*) as counts 
                             FROM %(track)s_binding  
                             WHERE CAST(pattern AS INT) != 0
                             GROUP BY pattern""" )


class BindingSummary(IntervalTracker):

    '''output summary counts of binding patterns.'''
    pattern = "(.*)_binding"

    def __call__(self, track):

        data = odict()

        data["binding"] = self.getValue("""
        SELECT COUNT(*) FROM %(track)s_binding
        WHERE overlap > 0""" )
        cols = [x for x in self.getColumns(
            "%s_binding" % track) if x.endswith("_overlap")]

        for section in ("flank5", "utr5", "cds", "first_exon", "first_intron", "middle_intron", "last_intron", "utr3", "flank3", "intron"):
            columns = [x for x in cols if x.startswith(section)]
            columns.sort()
            if section.endswith("5"):
                columns.reverse()

            for column in columns:
                data[column] = self.getValue( """
                              SELECT COUNT(*) FROM %(track)s_binding 
                              WHERE %(column)s > 0""" )
        return data


class BindingFullSummary(IntervalTracker):

    '''output summary counts of binding patterns.'''
    pattern = "(.*)_binding"

    def __call__(self, track):

        data = odict()

        # data["no binding"]  = self.getValue( """
        #                       SELECT COUNT(*) FROM %(track)s_binding
        #                       WHERE CAST(pattern AS INT) == 0""" )

        for section in ("flank5", "utr5", "cds", "first_intron", "middle_intron", "last_intron", "utr3", "flank3"):
            data[section] = self.getValue( """
                              SELECT COUNT(*) FROM %(track)s_binding 
                              WHERE %(section)s_overlap > 0""" )
        return data
