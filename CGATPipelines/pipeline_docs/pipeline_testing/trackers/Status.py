from TestingReport import *


class ComparisonStatus(Status):
    '''pipeline status
    '''

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM md5_compare")
        return tuple([x[0] for x in d])

    slices = ('MD5Comparison',)

    def testMD5Comparison(self, track):
        '''
        PASS: All files exist and are the same.

        FAIL: There are extra files or missing files

        WARN: All files exist, but there are differences in the files.

        The value indicates the number of files missing, exta
        or different.
        '''
        data = dict(self.getRow("""SELECT missing, extra, different
        FROM md5_compare WHERE track = '%(track)s'"""))

        total = sum(data.values())
        if total == 0:
            status = "PASS"
        elif data['missing'] == 0 and data['extra'] == 0:
            status = "WARN"
        else:
            status = "FAIL"

        value = ",".join(["%s:%i" % (x, y) for x, y in data.items()])

        return status, value
