from . import cnestedlist

import sqlite3
import os
import sys

if sys.version_info.major >= 3:
    import pickle as pickle
else:
    import cPickle as pickle

sqlite3.register_converter("pickle", pickle.loads)


class NCLSimple(object):
    """a nested contained list in memory storing
    no additional data.

    If *filestem* is not given, an new database is created
    in memory. If *filestem* is given, a database is opened 
    from disk. If it does not exist or *force* is set to True, 
    a new database is created. 

    Convenience wrapper around low-level ncl functions.

    Note that NCL objects are most efficient
    if used in write-once/read many applications.
    """

    def __init__(self, filestem=None, force=False):
        self.mTuples = []
        self.mIsDirty = False
        if filestem != None:
            self.mFilestem = filestem
            if not force and os.path.exists(
                    os.path.abspath(filestem) + ".idb"):
                self.mFromDisk = True
                self.mDatabase = cnestedlist.IntervalFileDB(filestem)
            else:
                self.mFromDisk = False
                self.mDatabase = cnestedlist.IntervalDB()
        else:
            self.mFilestem = None
            self.mFromDisk = False
            self.mDatabase = cnestedlist.IntervalDB()

    def add(self, start, end):
        """add segment *start*,*end* with value to database.

        returns the index of the added segment.
        """
        assert self.mFromDisk is False, "can not add to pre-existing or flushed databases"
        if start < 0:
            raise ValueError("only positive coordinates are accepted (%i<0)" % start)
        if start >= end:
            raise ValueError( "adding empty/invalid interval (%i,%i)" % (start,end))
        v = len(self.mTuples)
        self.mTuples.append((start, end, v))
        self.mIsDirty = True
        return v

    def find(self, start, end):
        """find intervals in database overlapping with *start* and *end*.

        returns an :class:`ncl.IntervalDBIterator`
        """
        if start < 0:
            raise ValueError("only positive coordinates are accepted (%i<0)" % start)
        self._commit()
        return self.mDatabase.find_overlap(start, end)

    def _commit(self):
        """commit database if changed."""
        if self.mIsDirty:
            self.mDatabase.fromlist(self.mTuples)
            self.mIsDirty = False

    def __del__(self):
        """flush database to disk."""
        if self.mFilestem and not self.mFromDisk:
            if self.mIsDirty:
                self.mDatabase.fromlist(self.mTuples)
            # flush database
            self.mDatabase.write_binaries(self.mFilestem)


class NCL(NCLSimple):
    """a nested contained list in memory. This
    class stores information on intervals.

    Note that NCL objects are most efficient
    if used in write-once/read many applications.
    """

    def __init__(self, *args, **kwargs):
        NCLSimple.__init__(self, *args, **kwargs)
        self.mValues = []
        # route calls to __getitem__ directly to list if in memory
        if self.mFromDisk:
            fn = self.mFilestem + ".vals" 
            self.mDBHandle = sqlite3.connect(fn)

    def add(self, start, end, value):
        """add segment *start*,*end* with value to database.

        returns the index of the added segment.
        """
        assert start >= 0, "only positive coordinates are accepted (%i<0)" % start
        if start >= end:
            raise ValueError(
                "adding empty/invalid interval (%i,%i)" % (start, end))
        self.mValues.append(value)
        return NCLSimple.add(self, start, end)

    def __getitem__(self, key):
        """get a value from the database
        """
        cc = self.mDBHandle.cursor()
        val = cc.execute(
            "SELECT value FROM data WHERE id = '%i'" % key).fetchone()[0]
        cc.close()
        return pickle.loads(str(val))

    def find(self, start, end):
        """find intervals overlapping *start* and *end*.

        returns an :class:`ncl.IteratorWithValues`
        """
        if self.mFromDisk:
            return IteratorWithValues(self, NCLSimple.find(self, start, end))
        else:
            return IteratorWithValues(self.mValues, NCLSimple.find(self, start, end))

    def _commit(self):
        """commit database if changed."""
        if self.mIsDirty and self.mFilestem:
            # flush database
            self._flushValues()

        NCLSimple._commit(self)

    def _flushValues(self):
        """flush values to disk."""

        fn = self.mFilestem + ".vals"
        if os.path.exists(fn):
            os.remove(fn)
        dbhandle = sqlite3.connect(fn)

        cc = dbhandle.cursor()
        cc.execute("""create table data( id INTEGER, value BLOB);""")
        cc.executemany(
            """INSERT INTO data VALUES (?,?)""",
            list(enumerate(map(lambda x: sqlite3.Binary(pickle.dumps(x)),
                               self.mValues))))
        dbhandle.commit()
        cc.close()

    def __del__(self):
        """flush database to disk."""
        if self.mFilestem and self.mIsDirty:
            self._flushValues()

        NCLSimple.__del__( self)


class IteratorWithValues(object):
    """an iterator over intervals with values.
    """

    def __init__(self, values, iter):

        self.mValues = values
        self.mIterator = iter

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        start, end, idx = self.mIterator.next()
        return (start, end, self.mValues[idx])
