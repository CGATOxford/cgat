'''IOTools - tools for I/O operations
==================================

This module contains utility functions for reading/writing from files.
These include methods for

* inspecting files, such as :func:`getFirstLine`, :func:`getLastLine`
  and :func:`isEmpty`,

* working with filenames, such as :func:`which` and :func:`snip`,
  :func:`checkPresenceOfFiles`

* manipulating file, such as :func:`openFile`, :func:`zapFile`,
  :func:`cloneFile`, :func:`touchFile`, :func:`shadowFile`.

* converting values for input/output, such as :func:`val2str`,
  :func:`str2val`, :func:`prettyPercent`, :func:`human2bytes`,
  :func:`convertDictionary`.

* iterating over file contents, such as :func:`iterate`,
  :func:`iterator_split`,

* creating lists/dictionaries from files, such as :func:`readMap` and
  :func:`readList`, and

* working with file collections (see :class:`FilePool`).

Reference
---------

'''

import collections
import glob
import gzip
import itertools
import numpy
import numpy.ma
import os
import re
import shutil
import stat
import string
import subprocess
import sys
import time


def getFirstLine(filename, nlines=1):
    """return the first line of a file.

    Arguments
    ---------
    filename : string
       The name of the file to be opened.
    nlines : int
       Number of lines to return.

    Returns
    -------
    string
       The first line(s) of the file.

    """
    # U is to open it with Universal newline support
    with open(filename, 'rU') as f:
        line = "".join([f.readline() for x in range(nlines)])
    return line


def getLastLine(filename, nlines=1, read_size=1024, encoding="utf-8"):
    """return the last line of a file.

    This method works by working back in blocks of `read_size` until
    the beginning of the last line is reached.

    Arguments
    ---------
    filename : string
       Name of the file to be opened.
    nlines : int
       Number of lines to return.
    read_size : int
       Number of bytes to read.

    Returns
    -------
    string
       The last line(s) of the file.

    """

    # py3 requires binary mode for negative seeks
    f = open(filename, 'rb')
    offset = read_size
    f.seek(0, 2)
    file_size = f.tell()
    if file_size == 0:
        return ""
    while 1:
        if file_size < offset:
            offset = file_size
        f.seek(-1 * offset, 2)
        read_str = f.read(offset)
        read_str = read_str.decode(encoding)
        lines = read_str.strip().splitlines()
        if len(lines) >= nlines + 1:
            return "\n".join(lines[-nlines:])
        if offset == file_size:   # reached the beginning
            return read_str
        offset += read_size
    f.close()


def getNumLines(filename, ignore_comments=True):
    """count number of lines in filename.

    Arguments
    ---------
    filename : string
       Name of the file to be opened.
    ignore_comments : bool
       If true, ignore lines starting with ``#``.

    Returns
    -------
    int
       The number of line(s) in the file.

    """

    if ignore_comments:
        filter_cmd = '| grep -v "#" '
    else:
        filter_cmd = ""

    # the implementation below seems to fastest
    # see https://gist.github.com/0ac760859e614cd03652
    # and
    # http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    if filename.endswith(".gz"):
        cmd = "zcat %(filename)s %(filter_cmd)s | wc -l" % locals()
    else:
        cmd = "cat %(filename)s %(filter_cmd)s | wc -l" % locals()

    out = subprocess.Popen(cmd,
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT
                           ).communicate()[0]
    return int(out.partition(b' ')[0])


def isEmpty(filename):
    """return True if file exists and is empty.

    Raises
    ------
    OSError
       If file does not exist
    """
    # don't now about stdin
    if filename == "-":
        return False
    return os.stat(filename)[stat.ST_SIZE] == 0


def isComplete(filename):
    '''return True if file exists and is complete.

    A file is complete if its last line starts
    with ``# job finished``.
    '''
    if filename.endswith(".gz"):
        raise NotImplementedError(
            'isComplete not implemented for compressed files')
    if isEmpty(filename):
        return False
    lastline = getLastLine(filename)
    return lastline.startswith("# job finished")


def touchFile(filename, times=None):
    '''update/create a sentinel file.

    Compressed files (ending in .gz) are created
    as empty 'gzip' files, i.e., with a header.
    '''
    existed = os.path.exists(filename)

    if filename.endswith(".gz") and not existed:
        # this will automatically add a gzip header
        fhandle = open(filename, 'a+b')
        fhandle = gzip.GzipFile(filename, fileobj=fhandle)
    else:
        fhandle = open(filename, 'a')

    try:
        os.utime(filename, times)
    finally:
        fhandle.close()


def openFile(filename, mode="r", create_dir=False, encoding="utf-8"):
    '''open file called *filename* with mode *mode*.

    gzip - compressed files are recognized by the
    suffix ``.gz`` and opened transparently.

    Note that there are differences in the file
    like objects returned, for example in the
    ability to seek.

    Arguments
    ---------
    filename : string
    mode : string
       File opening mode
    create_dir : bool
       If True, the directory containing filename
       will be created if it does not exist.

    Returns
    -------
    File or file-like object in case of gzip compressed files.
    '''

    _, ext = os.path.splitext(filename)

    if create_dir:
        dirname = os.path.dirname(filename)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)

    if ext.lower() in (".gz", ".z"):
        if sys.version_info.major >= 3:
            if mode == "r":
                return gzip.open(filename, 'rt', encoding=encoding)
            elif mode == "w":
                return gzip.open(filename, 'wt', encoding=encoding)
            elif mode == "a":
                return gzip.open(filename, 'wt', encoding=encoding)
        else:
            return gzip.open(filename, mode)
    else:
        return open(filename, mode, encoding=encoding)


def force_str(iterator, encoding="ascii"):
    """iterate over lines in iterator and force to string"""
    if sys.version_info.major >= 3:
        for line in iterator:
            yield line.decode(encoding)
    else:
        for line in iterator:
            yield line


def zapFile(filename, outfile=None):
    '''replace *filename* with empty file.

    File attributes such as accession times are preserved.

    If the file is a link, the link will be broken and replaced with
    an empty file having the same attributes as the file linked to.

    It also takes an optional outfile. If the outfile has zero byte,
        it usually means there's an error in generating the outfile,
        and it will throw an error and stop.

    Returns
    -------
    stat_object
       A stat object of the file cleaned.
    link_destination : string
       If the file was a link, the file being linked to.

    '''
    # outfile as zero byte? Let's throw an error and stop
    if outfile and os.path.getsize(outfile) == 0:
        raise ValueError('%s has size zero!' % outfile)

    # stat follows times to links
    original = os.stat(filename)

    # return if file already has size 0
    if original.st_size == 0:
        return None, None

    if os.path.islink(filename):
        linkdest = os.readlink(filename)
        os.unlink(filename)
        f = open(filename, "w")
        f.close()
    else:
        linkdest = None
        f = open(filename, "w")
        f.truncate()
        f.close()

    # Set original times
    os.utime(filename, (original.st_atime, original.st_mtime))
    os.chmod(filename, original.st_mode)

    return original, linkdest


def cloneFile(infile, outfile):
    '''create a clone of ``infile`` named ``outfile``
    by creating a soft-link.
    '''
    # link via relative paths, otherwise it
    # fails if infile and outfile are in different
    # directories or in a subdirectory
    if os.path.dirname(infile) != os.path.dirname(outfile):
        relpath = os.path.relpath(
            os.path.dirname(infile), os.path.dirname(outfile))
    else:
        relpath = "."
    target = os.path.join(relpath, os.path.basename(infile))

    try:
        os.symlink(target, outfile)
    except OSError:
        pass


def shadowFile(infile, outfile):
    '''move ```infile``` as ```outfile```, and
    touch ```infile```.
    This could be useful when one wants to skip
    some steps in a pipeline.
    Note that zapFile is not needed when shadowFile
    is used
    '''
    if outfile != infile:
        shutil.move(infile, outfile)
        touchFile(infile)
        # reset outfile's timestamp
        time.sleep(1)
        touchFile(outfile)
    else:
        raise ValueError('Panic: infile and outfile names cannot be the same')


def val2str(val, format="%5.2f", na="na"):
    '''return a formatted value.

    If value does not fit format string, return "na"
    '''
    if type(val) == int:
        return format % val
    elif type(val) == float:
        return format % val

    try:
        x = format % val
    except (ValueError, TypeError):
        x = na
    return x


def str2val(val, format="%5.2f", na="na", list_detection=False):
    """guess type (int, float) of value.

    If `val` is neither int nor float, the value
    itself is returned.
    """

    if val is None:
        return val

    def _convert(v):
        try:
            x = int(v)
        except ValueError:
            try:
                x = float(v)
            except ValueError:
                return v
        return x

    if list_detection and "," in val:
        return [_convert(v) for v in val.split(",")]
    else:
        return _convert(val)


def prettyPercent(numerator, denominator, format="%5.2f", na="na"):
    """output a percent value or "na" if not defined"""
    try:
        x = format % (100.0 * numerator / denominator)
    except (ValueError, ZeroDivisionError):
        x = "na"
    return x


def prettyString(val):
    '''output val or na if val is None'''
    if val is not None:
        return val
    else:
        return "na"


def which(program):
    """check if `program` is in PATH and is executable.

    Returns
    -------
    string
       The full path to the program. Returns None if not found.

    """
    # see http://stackoverflow.com/questions/377017/test-if-
    #  executable-exists-in-python

    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def iterate(infile):
    '''iterate over infile and return a :py:class:`collections.namedtuple`
    according to a header in the first row.

    Lines starting with ``#`` are skipped.

    '''

    n = 0
    for line in infile:
        if line.startswith("#"):
            continue
        n += 1
        if n == 1:
            # replace non-alphanumeric characters with _
            header = re.sub("[^a-zA-Z0-9_\s]", "_", line[:-1]).split()
            DATA = collections.namedtuple("DATA", header)
            continue

        result = DATA(*line[:-1].split())

        yield result


def iterate_tabular(infile, sep="\t"):
    '''iterate over file `infile` skipping lines starting with
    ``#``.

    Within a line, records are separated by `sep`.

    Yields
    ------
    tuple
        Records within a line

    '''
    for line in infile:
        if line.startswith("#"):
            continue
        yield line[:-1].split(sep)


def iterator_split(infile, regex):
    '''Return an iterator of file chunks based on a known logical start
    point `regex` that splits the file into intuitive chunks.  This
    assumes the file is structured in some fashion.  For arbitrary
    number of bytes use file.read(`bytes`).  If a header is present it
    is returned as the first file chunk.

    infile must be either an open file handle or an iterable.

    '''
    chunk_list = []

    regex = re.compile(regex)

    for x in infile:
        if regex.search(x):
            if len(chunk_list):
                # return the current chunk and start a new one from this point
                yield chunk_list
            chunk_list = []
            chunk_list.append(x)
        else:
            chunk_list.append(x)
    yield chunk_list


def snip(filename, extension=None, alt_extension=None,
         strip_path=False):
    '''return prefix of `filename`, that is the part without the
    extension.

    If `extension` is given, make sure that filename has the
    extension (or `alt_extension`).

    If `strip_path` is set to true, the path is stripped from the file
    name.

    '''
    if extension:
        if filename.endswith(extension):
            root = filename[:-len(extension)]
        elif alt_extension and filename.endswith(alt_extension):
            root = filename[:-len(alt_extension)]
        else:
            raise ValueError("'%s' expected to end in '%s'" %
                             (filename, extension))
    else:
        root, ext = os.path.splitext(filename)

    if strip_path:
        snipped = os.path.basename(root)
    else:
        snipped = root

    return snipped


def checkPresenceOfFiles(filenames):
    """check for the presence/absence of files

    Parameters
    ----------
    filenames : list
        Filenames to check for presence.

    Returns
    -------
    missing : list
        List of missing filenames
    """

    missing = []
    for filename in filenames:
        if not os.path.exists(filename):
            missing.append(filename)
    return missing


def human2bytes(s):
    """
    Attempts to guess the string format based on default symbols
    set and return the corresponding bytes as an integer.
    When unable to recognize the format ValueError is raised.

      >>> human2bytes('0 B')
      0
      >>> human2bytes('1 K')
      1024
      >>> human2bytes('1 M')
      1048576
      >>> human2bytes('1 Gi')
      1073741824
      >>> human2bytes('1 tera')
      1099511627776

      >>> human2bytes('0.5kilo')
      512
      >>> human2bytes('0.1  byte')
      0
      >>> human2bytes('1 k')  # k is an alias for K
      1024
      >>> human2bytes('12 foo')
      Traceback (most recent call last):
          ...
      ValueError: can't interpret '12 foo'

    Author: Giampaolo Rodola' <g.rodola [AT] gmail [DOT] com>
    License: MIT

    https://gist.github.com/leepro/9694638
    """
    SYMBOLS = {
        'customary': ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
        'customary_ext': ('byte', 'kilo', 'mega', 'giga',
                          'tera', 'peta', 'exa',
                          'zetta', 'iotta'),
        'iec': ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
        'iec_ext': ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi',
                    'zebi', 'yobi'),
    }

    init = s
    num = ""
    while s and s[0:1].isdigit() or s[0:1] == '.':
        num += s[0]
        s = s[1:]
    num = float(num)
    letter = s.strip()
    for name, sset in list(SYMBOLS.items()):
        if letter in sset:
            break
    else:
        if letter == 'k':
            # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
            sset = SYMBOLS['customary']
            letter = letter.upper()
        else:
            raise ValueError("can't interpret %r" % init)
    prefix = {sset[0]: 1}
    for i, s in enumerate(sset[1:]):
        prefix[s] = 1 << (i + 1) * 10

    return int(num * prefix[letter])


def convertDictionary(d, map={}):
    """convert string values in a dictionary to numeric types.

    Arguments
    d : dict
       The dictionary to convert
    map : dict
       If map contains 'default', a default conversion is enforced.
       For example, to force int for every column but column ``id``,
       supply map = {'default' : "int", "id" : "str" }
    """

    rx_int = re.compile("^\s*[+-]*[0-9]+\s*$")
    rx_float = re.compile("^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")

    # pre-process with 'default'
    if "default" in map:
        k = "default"
        if map[k] == "int":
            default = int
        elif map[k] == "float":
            default = float
        elif map[k] == "string":
            default = str
    else:
        default = False

    for k, vv in list(d.items()):

        if vv is None:
            continue
        v = vv.strip()
        try:
            if k in map:
                if map[k] == "int":
                    d[k] = int(v)
                elif map[k] == "float":
                    d[k] = float(v)
                elif map[k] == "string":
                    pass
                continue
            elif default:
                if v != "":
                    d[k] = default(v)
                else:
                    d[k] = v
                continue
        except TypeError as msg:
            raise TypeError("conversion in field: %s, %s" % (k, msg))

        try:
            if rx_int.match(v):
                d[k] = int(v)
            elif rx_float.match(v):
                d[k] = float(v)
        except TypeError as msg:
            raise TypeError(
                "expected string or buffer: offending value = '%s' " % str(v))
        except ValueError as msg:
            raise ValueError("conversion error: %s, %s" % (msg, str(d)))
    return d


class nested_dict(collections.defaultdict):
    """Auto-vivifying nested dictionaries.

    For example::

      nd= nested_dict()
      nd["mouse"]["chr1"]["+"] = 311

   """

    def __init__(self):
        collections.defaultdict.__init__(self, nested_dict)

    def iterflattened(self):
        """
        iterate through values with nested keys flattened into a tuple
        """

        for key, value in self.items():
            if isinstance(value, nested_dict):
                for keykey, value in value.iterflattened():
                    yield (key,) + keykey, value
            else:
                yield (key,), value


def flatten(l, ltypes=(list, tuple)):
    '''flatten a nested list.

    This method works with any list-like container
    such as tuples.

    Arguments
    ---------
    l : list
        A nested list.
    ltypes : list
        A list of valid container types.

    Returns
    -------
    list : list
        A flattened list.
    '''
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def invert_dictionary(dict, make_unique=False):
    """returns an inverted dictionary with keys and values swapped.
    """
    inv = {}
    if make_unique:
        for k, v in dict.items():
            inv[v] = k
    else:
        for k, v in dict.items():
            inv.setdefault(v, []).append(k)
    return inv


class FilePool:
    """manage a pool of output files.

    This class will keep a large number of files open. To
    see if you can handle this, check the limit within the shell::

       ulimit -n

    The number of currently open and maximum open files in the system:

      cat /proc/sys/fs/file-nr

    Changing these limits might not be easy without root privileges.

    The maximum number of files opened is given by :attr:`maxopen`.
    This class is inefficient if the number of files is larger than
    :attr:`maxopen` and calls to `write` do not group keys together.

    To use this class, create a FilePool and write to it as if it was
    a single file, specifying a section for each write::

        pool = FilePool("%s.tsv")
        for value in range(100):
            for section in ("file1", "file2", "file3"):
                 pool.write(section, str(value) + ",")

    This will create three files called ``file1.tsv``, ``file2.tsv``,
    ``file3.tsv``, each containing the numbers from 0 to 99.

    The FilePool acts otherwise as a dictionary providing access to
    the number of times an item has been written to each file::

        print pool["file1]
        print pool.items()

    Parameters
    ----------

    output_pattern : string
       output pattern to use. Should contain a "%s". If set to None, the
       pattern "%s" will be used.
    header : string
       optional header to write when writing to a file the first time.
    force : bool
       overwrite existing files. All files matching the pattern will be
       deleted.

    """

    maxopen = 5000

    def __init__(self,
                 output_pattern=None,
                 header=None,
                 force=True):

        self.mFiles = {}
        self.mOutputPattern = output_pattern
        self.mCounts = collections.defaultdict(int)
        self.mHeader = header
        if force and output_pattern:
            for f in glob.glob(re.sub("%s", "*", output_pattern)):
                os.remove(f)

    def __del__(self):
        """close all open files."""
        for file in list(self.mFiles.values()):
            file.close()

    def __len__(self):
        return len(self.mCounts)

    def close(self):
        """close all open files."""
        for file in list(self.mFiles.values()):
            file.close()

    def values(self):
        return list(self.mCounts.values())

    def keys(self):
        return list(self.mCounts.keys())

    def iteritems(self):
        return iter(self.mCounts.items())

    def items(self):
        return list(self.mCounts.items())

    def __iter__(self):
        return self.mCounts.__iter__()

    def getFile(self, identifier):
        return identifier

    def getFilename(self, identifier):
        """get filename for an identifier."""

        if self.mOutputPattern:
            return re.sub("%s", str(identifier), self.mOutputPattern)
        else:
            return identifier

    def setHeader(self, header):
        """set the header to be written to each file when opening
        for the first time."""

        self.mHeader = header

    def openFile(self, filename, mode="w"):
        """open file.

        If file is in a new directory, create directories.
        """
        if mode in ("w", "a"):
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)

        return openFile(filename, mode)

    def write(self, identifier, line):
        """write `line` to file specified by `identifier`"""
        filename = self.getFilename(identifier)

        if filename not in self.mFiles:

            if self.maxopen and len(self.mFiles) > self.maxopen:
                for f in list(self.mFiles.values()):
                    f.close()
                self.mFiles = {}

            self.mFiles[filename] = openFile(filename, "a")
            if self.mHeader:
                self.mFiles[filename].write(self.mHeader)

        try:
            self.mFiles[filename].write(line)
        except ValueError as msg:
            raise ValueError(
                "error while writing to %s: msg=%s" % (filename, msg))
        self.mCounts[filename] += 1

    def deleteFiles(self, min_size=0):
        """delete all files below a minimum size `min_size` bytes."""

        ndeleted = 0
        for filename, counts in list(self.mCounts.items()):
            if counts < min_size:
                os.remove(filename)
                ndeleted += 1

        return ndeleted


class FilePoolMemory(FilePool):
    """manage a pool of output files in memory.

    The usage is the same as :class:`FilePool` but the data is cached
    in memory before writing to disk.

    """

    maxopen = 5000

    def __init__(self, *args, **kwargs):
        FilePool.__init__(self, *args, **kwargs)

        self.data = collections.defaultdict(list)
        self.isClosed = False

    def __del__(self):
        """close all open files.
        """
        if not self.isClosed:
            self.close()

    def close(self):
        """close all open files.
        writes the data to disk.
        """
        if self.isClosed:
            raise IOError("write on closed FilePool in close()")

        for filename, data in self.data.items():
            f = openFile(filename, "a")
            if self.mHeader:
                f.write(self.mHeader)
            f.write("".join(data))
            f.close()

        self.isClosed = True

    def write(self, identifier, line):

        filename = self.getFilename(identifier)
        self.data[filename].append(line)
        self.mCounts[filename] += 1


def readMap(infile,
            columns=(0, 1),
            map_functions=(str, str),
            both_directions=False,
            has_header=False,
            dtype=dict):
    """read a map (key, value pairs) from infile.

    If there are multiple entries for the same key, only the
    last entry will be recorded.

    Arguments
    ---------
    infile : File
       File object to read from
    columns : tuple
       Columns (A, B) to take from the file to create the mapping from
       A to B.
    map_functions : tuple
       Functions to convert the values in the rows to the desired
       object types such as int or float.
    both_directions : bool
       If true, both mapping directions are returned.
    has_header : bool
       If true, ignore first line with header.
    dtype : function
       datatype to use for the dictionaries.

    Returns
    -------
    map : dict
       A dictionary containing the mapping. If `both_directions` is true,
       two dictionaries will be returned.

    """
    m = dtype()
    r = dtype()
    n = 0

    if columns == "all":
        key_column = 0
        value_column = None
    else:
        key_column, value_column = columns

    key_function, value_function = map_functions
    # default is to return a tuple for multiple values
    datatype = None

    for l in infile:
        if l[0] == "#":
            continue
        n += 1

        if has_header and n == 1:
            if columns == "all":
                header = l[:-1].split("\t")
                # remove the first column
                datatype = collections.namedtuple("DATA", header[1:])
            continue

        d = l[:-1].split("\t")
        if len(d) < 2:
            continue
        key = key_function(d[key_column])
        if value_column:
            val = value_function(d[value_column])
        elif datatype:
            val = datatype._make([d[x] for x in range(1, len(d))])
        else:
            val = tuple(map(value_function, [d[x] for x in range(1, len(d))]))

        m[key] = val
        if val not in r:
            r[val] = []
        r[val].append(key)

    if both_directions:
        return m, r
    else:
        return m


def readList(infile,
             column=0,
             map_function=str,
             map_category={},
             with_title=False):
    """read a list of values from infile.

    Arguments
    ---------
    infile : File
       File object to read from
    columns : int
       Column to take from the file.
    map_function : function
       Function to convert the values in the rows to the desired
       object types such as int or float.
    map_category : dict
       When given, automatically transform/map the values given
       this dictionary.
    with_title : bool
       If true, first line of file is title and will be ignored.

    Returns
    -------
    list : list
       A list with the values.
    """

    m = []
    title = None
    for l in infile:
        if l[0] == "#":
            continue
        if with_title and not title:
            title = l[:-1].split("\t")[column]
            continue

        try:
            d = map_function(l[:-1].split("\t")[column])
        except ValueError:
            continue

        if map_category:
            d = map_category[d]
        m.append(d)

    return m


def readMultiMap(infile,
                 columns=(0, 1),
                 map_functions=(str, str),
                 both_directions=False,
                 has_header=False,
                 dtype=dict):
    """read a map (pairs of values) from infile.

    In contrast to :func:`readMap`, this method permits multiple
    entries for the same key.

    Arguments
    ---------
    infile : File
       File object to read from
    columns : tuple
       Columns (A, B) to take from the file to create the mapping from
       A to B.
    map_functions : tuple
       Functions to convert the values in the rows to the desired
       object types such as int or float.
    both_directions : bool
       If true, both mapping directions are returned in a tuple, i.e.,
       A->B and B->A.
    has_header : bool
       If true, ignore first line with header.
    dtype : function
       datatype to use for the dictionaries.

    Returns
    -------
    map : dict
       A dictionary containing the mapping. If `both_directions` is true,
       two dictionaries will be returned.

    """
    m = dtype()
    r = dtype()
    n = 0
    for l in infile:
        if l[0] == "#":
            continue
        n += 1

        if has_header and n == 1:
            continue

        d = l[:-1].split("\t")
        try:
            key = map_functions[0](d[columns[0]])
            val = map_functions[1](d[columns[1]])
        except (ValueError, IndexError) as msg:
            raise ValueError("parsing error in line %s: %s" % (l[:-1], msg))

        if key not in m:
            m[key] = []
        m[key].append(val)
        if val not in r:
            r[val] = []
        r[val].append(key)

    if both_directions:
        return m, r
    else:
        return m


def readMatrix(infile, dtype=numpy.float):
    '''read a numpy matrix from infile.

    return tuple of matrix, row_headers, col_headers
    '''

    lines = [l for l in infile.readlines() if not l.startswith("#")]
    nrows = len(lines) - 1
    col_headers = lines[0][:-1].split("\t")[1:]
    ncols = len(col_headers)
    matrix = numpy.zeros((nrows, ncols), dtype=dtype)
    row_headers = []

    for row, l in enumerate(lines[1:]):
        data = l[:-1].split("\t")
        row_headers.append(data[0])
        matrix[row] = numpy.array(data[1:], dtype=dtype)

    return matrix, row_headers, col_headers


def writeMatrix(outfile, matrix, row_headers, col_headers, row_header=""):
    '''write a numpy matrix to outfile.

    *row_header* gives the title of the rows
    '''

    outfile.write("%s\t%s\n" % (row_header, "\t".join(col_headers)))
    for x, row in enumerate(matrix):
        assert len(row) == len(col_headers)
        outfile.write("%s\t%s\n" % (row_headers[x], "\t".join(map(str, row))))


def readTable(file,
              separator="\t",
              numeric_type=numpy.float,
              take="all",
              headers=True,
              truncate=None,
              cumulate_out_of_range=True,
              ):
    """read a table of values.

    If cumulate_out_of_range is set to true, the terminal bins will
    contain the cumulative values of bins out of range.

    .. note:: Deprecated
       use pandas dataframes instead

    """

    lines = [x for x in file.readlines() if x[0] != "#"]

    if len(lines) == 0:
        return None, []

    if take == "all":
        num_cols = len(string.split(lines[0][:-1], "\t"))
        take = list(range(0, num_cols))
    else:
        num_cols = len(take)

    if headers:
        headers = lines[0][:-1].split("\t")
        headers = [headers[x] for x in take]
        del lines[0]

    num_rows = len(lines)
    matrix = numpy.ma.masked_array(
        numpy.zeros((num_rows, num_cols), numeric_type))

    if truncate:
        min_row, max_row = truncate

    nrow = 0
    min_data = [0] * num_cols
    max_data = None
    for l in lines:
        data = l[:-1].split("\t")
        data = [data[x] for x in take]

        # try conversion. Unparseable fields set to missing_value
        for x in range(len(data)):
            try:
                data[x] = float(data[x])
            except ValueError:
                data[x] = numpy.ma.masked

        if truncate is not None:
            if data[0] < min_row:
                if cumulate_out_of_range:
                    for x in range(1, num_cols):
                        min_data[x] += data[x]
                continue
            elif data[0] >= max_row:
                if max_data is None:
                    max_data = [0] * num_cols
                    max_data[0] = max_row
                for x in range(1, num_cols):
                    try:
                        max_data[x] += data[x]
                    except TypeError:
                        # missing values cause type errors
                        continue
                continue
            elif min_row is not None:
                if cumulate_out_of_range:
                    for x in range(0, num_cols):
                        try:
                            min_data[x] += data[x]
                        except TypeError:
                            # missing values cause type errors
                            continue
                    else:
                        min_data = data
                data = min_data
                min_row = None

        # copy values into matrix
        # this is a bit clumsy, but missing values
        # cause an error otherwise
        for x in range(len(data)):
            matrix[nrow, x] = data[x]

        nrow += 1

    if truncate is not None:
        if cumulate_out_of_range:
            if max_data is not None:
                matrix[nrow] = max_data

        # truncate matrix
        matrix = matrix[0:nrow + 1, 0:num_cols]

    return matrix, headers


def writeTable(outfile, table, columns=None, fillvalue=""):
    '''write a table to outfile.

    If table is a dictionary, output columnwise. If *columns* is a list,
    only output columns in columns in the specified order.

    .. note:: Deprecated
       use pandas dataframes instead

    '''

    if type(table) == dict:
        if columns is None:
            columns = list(table.keys())
        outfile.write("\t".join(columns) + "\n")
        # get data
        data = [table[x] for x in columns]
        # transpose
        data = list(itertools.zip_longest(*data, fillvalue=fillvalue))

        for d in data:
            outfile.write("\t".join(map(str, d)) + "\n")

    else:
        raise NotImplementedError


def ReadMap(*args, **kwargs):
    """deprecated, use readMap."""
    return readMap(*args, **kwargs)


def ReadList(*args, **kwargs):
    """deprecated, use readList()"""
    return readList(*args, **kwargs)


def writeLines(outfile, lines, header=False):
    ''' expects [[[line1-field1],[line1-field2 ] ],... ]'''
    handle = openFile(outfile, "w")

    if header:
        handle.write("\t".join([str(title) for title in header]) + "\n")

    for line in lines:
        handle.write("\t".join([str(field) for field in line]) + "\n")

    handle.close()


def txtToDict(filename, key=None, sep="\t"):
    '''make a dictionary from a text file keyed
    on the specified column.'''

    # Please see function in readDict()
    count = 0
    result = {}
    valueidx, keyidx = False, False
    field_names = []

    with open(filename, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if count == 0:
                fieldn = 0
                for rawfield in line.split(sep):
                    field = rawfield.strip()
                    if field == key:
                        keyidx = fieldn
                    field_names.append(field)
                    fieldn += 1

                if not keyidx:
                    raise ValueError("key name not found in header")
                # if not valueidx:
                #   raise ValueError(
                #     "value name not found in header")
            else:
                fields = [x.strip() for x in line.split(sep)]
                fieldn = 0
                thiskey = fields[keyidx]
                result[thiskey] = {}
                for field in fields:
                    if fieldn == keyidx:
                        pass
                    else:
                        colkey = field_names[fieldn]
                        result[thiskey][colkey] = field
                    fieldn += 1
            count += 1

    return(result)


def pickle(file_name, obj):
    '''dump a python object to a file using pickle'''
    with open(file_name, "wb") as pkl_file:
        pickle.dump(obj, pkl_file)
    return


def unpickle(file_name):
    '''retrieve a pickled python object from a file'''
    with open(file_name, "r") as pkl_file:
        data = pickle.load(pkl_file)
    return data
