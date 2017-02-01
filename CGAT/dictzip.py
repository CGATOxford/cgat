"""Functions that read and write gzipped files.

The user of the file doesn't have to worry about the compression.
Seeks are allowed, but are efficient only for files compressed by
the dictzip utility (which adhere to the gzip format).  Files may
be concatenated to overcome dictzip's 1.8 Gb size limit."""


import struct
import sys
import time
import bisect
import zlib
import builtins

__all__ = ["GzipFile", "open"]

FTEXT, FHCRC, FEXTRA, FNAME, FCOMMENT = 1, 2, 4, 8, 16

READ, WRITE = 1, 2

SEEK_SET, SEEK_CUR, SEEK_END = 0, 1, 2


def U32(i):
    """Return i as an unsigned integer, assuming it fits in 32 bits.

    If it's >= 2GB when viewed as a 32-bit unsigned int, return a long.
    """
    if i < 0:
        i += 1 << 32
    return i


def LOWU32(i):
    """Return the low-order 32 bits of an int, as a non-negative int."""
    return i & 0xFFFFFFFF


def write32(output, value):
    output.write(struct.pack("<l", value))


def write32u(output, value):
    # The L format writes the bit pattern correctly whether signed
    # or unsigned.
    output.write(struct.pack("<L", value))


def read32(input):
    return struct.unpack("<l", input.read(4))[0]


def open(filename, mode="rb", compresslevel=9, buffersize=None, chunksize=58315):
    """Shorthand for GzipFile(filename, mode, compresslevel, buffersize, chunksize).

    The filename argument is required; mode defaults to 'rb', compresslevel
    defaults to 9, buffersize to None (no random access points) and chunksize to 58315
    (good compression but slowish seeks; irrelevant without random access points)

    """
    return GzipFile(filename, mode, compresslevel, buffersize=buffersize, chunksize=chunksize)


class GzipFile:

    """The GzipFile class simulates most of the methods of a file object with
    the exception of the readinto() and truncate() methods.

    """

    myfileobj = None

    def __init__(self, filename=None, mode=None,
                 compresslevel=9, fileobj=None, buffersize=None, chunksize=58315):
        """Constructor for the GzipFile class.

        At least one of fileobj and filename must be given a
        non-trivial value.

        The new class instance is based on fileobj, which can be a regular
        file, a StringIO object, or any other object which simulates a file.
        It defaults to None, in which case filename is opened to provide
        a file object.

        When fileobj is not None, the filename argument is only used to be
        included in the gzip file header, which may includes the original
        filename of the uncompressed file.  It defaults to the filename of
        fileobj, if discernible; otherwise, it defaults to the empty string,
        and in this case the original filename is not included in the header.

        The mode argument can be any of 'r', 'rb', 'a', 'ab', 'w', or 'wb',
        depending on whether the file will be read or written.  The default
        is the mode of fileobj if discernible; otherwise, the default is 'rb'.
        Be aware that only the 'rb', 'ab', and 'wb' values should be used
        for cross-platform portability.

        The compresslevel argument is an integer from 1 to 9 controlling the
        level of compression; 1 is fastest and produces the least compression,
        and 9 is slowest and produces the most compression.  The default is 9.

        A nonzero buffersize argument instructs GZip to do buffered compression,
        allowing it to include a dictzip field in the header with flush points
        for random access.  The chunksize argument determines the distance between
        flush points; smaller values means faster random access but lower
        compression.  The default value is close to maximum compression.

        """

        # guarantee the file is opened in binary mode on platforms
        # that care about that sort of thing
        if mode and 'b' not in mode:
            mode += 'b'
        if fileobj is None:
            fileobj = self.myfileobj = builtins.open(filename, mode or 'rb')
        if filename is None:
            if hasattr(fileobj, 'name'):
                filename = fileobj.name
            else:
                filename = ''
        if mode is None:
            if hasattr(fileobj, 'mode'):
                mode = fileobj.mode
            else:
                mode = 'rb'

        if mode[0:1] == 'r':
            self.mode = READ
            # Set flag indicating start of a new member
            self._new_member = True
            # Set flag indicating normal gzip format
            self.dictzip = False
            self.extrabuf = ""
            self.extrasize = 0
            self.filename = filename

        elif mode[0:1] == 'w' or mode[0:1] == 'a':
            self.mode = WRITE
            if buffersize:
                self.dictzip = True
            else:
                self.dictzip = False
            self.compresslevel = compresslevel
            self.chunksize = chunksize
            self.buffersize = buffersize
            # dictzip's default chunk size of 58315 is too conservative
            if chunksize > 65400:
                raise IOError("Chunk size " +
                              str(chunksize) + " is too large; maximum is 65400")
            if self.dictzip and buffersize // chunksize > 32764:
                raise IOError("Buffer size " +
                              str(buffersize) +
                              " is too large; may result in too many chunks")

            self._init_write(filename)

        else:
            raise IOError("Mode " + mode + " not supported")

        self.fileobj = fileobj
        self.offset = 0

        if self.mode == WRITE:
            if self.dictzip:
                # intialize write buffer
                self.writebuf = ''
            else:
                # for ordinary gzip files, write header now
                self._write_gzip_header()

        if self.mode == READ:
            # read the headers of all dictzip members, to build database of
            # flush points
            # offset of member within uncompressed stream
            self.memberoffset = []
            self.memberchlen = []         # chunk length
            # absolute flush points within this member
            self.memberflushpoints = []
            self.dictzip = True
            pos = self.fileobj.tell()
            try:
                offset = 0
                while True:
                    self._iseof()
                    dictzipdata = self._read_gzip_header()
                    if dictzipdata:
                        chlen, flushpoints = dictzipdata
                        self.memberoffset.append(offset)
                        self.memberchlen.append(chlen)
                        for idx in range(len(flushpoints)):
                            flushpoints[idx] += self.fileobj.tell()
                        # keep flushpoints, including the one pointing beyond
                        # the data stream
                        self.memberflushpoints.append(flushpoints)
                        # point to length field at end of this member
                        # Add 4 to skip over the CRC32 field
                        # (I don't understand the "+2" -- header bytes?)
                        newpos = flushpoints[-1] + 2 + 4
                        self.fileobj.seek(newpos)
                        # will not exceed 2 Gb for dictzip files
                        isize = U32(read32(self.fileobj))
                        offset += isize
                    else:
                        self.dictzip = False
                        break
            except EOFError:
                pass
            self.uncompressed_length = offset
            self.fileobj.seek(pos)

    def __repr__(self):
        s = repr(self.fileobj)
        return '<gzip ' + s[1:-1] + ' ' + hex(id(self)) + '>'

    def _init_write(self, filename):
        if filename is not None and filename[-3:] != '.gz' and filename[-3:] != '.dz':
            if self.dictzip:
                filename = filename + '.dz'
            else:
                filename = filename + '.gz'
        self.filename = filename
        self.crc = zlib.crc32("")
        self.size = 0
        self.compress = zlib.compressobj(self.compresslevel,
                                         zlib.DEFLATED,
                                         -zlib.MAX_WBITS,
                                         zlib.DEF_MEM_LEVEL,
                                         0)

    def _write_gzip_header(self, size=None):
        self.fileobj.write('\037\213')             # magic header
        self.fileobj.write('\010')                 # compression method
        flags = 0
        if self.filename:
            flags = FNAME
        if self.dictzip:
            flags |= FEXTRA
        self.fileobj.write(chr(flags))
        write32u(self.fileobj, int(time.time()))
        self.fileobj.write('\002')                 # extraflag
        self.fileobj.write('\377')                 # os (unknown)
        if self.dictzip:
            chunks = 1 + (size - 1) // self.chunksize
            xlen = 10 + 2 * chunks
            # length of extra field
            self.fileobj.write(chr(xlen % 256))
            self.fileobj.write(chr(xlen // 256))
            # dictzip's magic word - 'R'andom 'A'ccess
            self.fileobj.write('RA')
            sublen = xlen - 4
            # length of subfield
            self.fileobj.write(chr(sublen % 256))
            self.fileobj.write(chr(sublen // 256))
            # dictzip header version
            self.fileobj.write('\001\000')
            self.fileobj.write(chr(self.chunksize % 256))    # size of chunk
            self.fileobj.write(chr(self.chunksize // 256))
            # number of chunks
            self.fileobj.write(chr(chunks % 256))
            self.fileobj.write(chr(chunks // 256))
            self.chunktablepos = self.fileobj.tell()
            for chunk in range(chunks):
                self.fileobj.write('\000\000')                  # placeholders
        if self.filename:
            self.fileobj.write(self.filename[:-3] + '\000')

    def _init_read(self):
        self.crc = zlib.crc32("")
        self.size = 0

    def _read_gzip_extra(self):
        xlen = ord(self.fileobj.read(1))
        xlen = xlen + 256 * ord(self.fileobj.read(1))
        xtra = self.fileobj.read(xlen)
        xptr = 0
        # loop over subfields
        while xptr < xlen:
            if xlen - xptr < 4:
                # ill-formed header: magic word + subfield length required
                return
            # subfield length
            sublen = ord(xtra[xptr + 2]) + 256 * ord(xtra[xptr + 3])
            ptr = xptr
            xptr += sublen + 4
            if xtra[ptr] != 'R' or xtra[ptr + 1] != 'A':
                continue     # magic word for dictzip data is 'R'andom 'A'ccess
            if xtra[ptr + 4] != '\001' or xtra[ptr + 5] != '\000':
                raise IOError("Unrecognized DictZip version: " +
                              str(ord(xtra[ptr + 4]) + 256 * ord(xtra[ptr + 5])))
            # chunk length
            chlen = ord(xtra[ptr + 6]) + 256 * ord(xtra[ptr + 7])
            # chunk count
            chcnt = ord(xtra[ptr + 8]) + 256 * ord(xtra[ptr + 9])
            if chcnt * 2 != sublen - 6:
                raise IOError("Invalid DictZip header: wrong number of chunks:" +
                              str(chcnt) + " expected " + str((sublen - 6) // 2))
            flushpoints = [0]
            for idx in range(chcnt):
                flushpoints.append(
                    flushpoints[-1] + ord(xtra[ptr + 10 + 2 * idx]) + 256 * ord(xtra[ptr + 11 + 2 * idx]))
            # ignore other subfields
            return (chlen, flushpoints)
        if xptr != xlen:
            raise IOError("Bad extra field in gzip header")
        return None

    def _read_gzip_header(self):
        magic = self.fileobj.read(2)
        if magic != '\037\213':
            raise IOError('Not a gzipped file')
        method = ord(self.fileobj.read(1))
        if method != 8:
            raise IOError('Unknown compression method')
        flag = ord(self.fileobj.read(1))
        # modtime = self.fileobj.read(4)
        # extraflag = self.fileobj.read(1)
        # os = self.fileobj.read(1)
        self.fileobj.read(6)
        dictzipdata = None

        if flag & FEXTRA:
            # Read the extra field
            dictzipdata = self._read_gzip_extra()
        if flag & FNAME:
            # Read and discard a null-terminated string containing the filename
            while True:
                s = self.fileobj.read(1)
                if not s or s == '\000':
                    break
        if flag & FCOMMENT:
            # Read and discard a null-terminated string containing a comment
            while True:
                s = self.fileobj.read(1)
                if not s or s == '\000':
                    break
        if flag & FHCRC:
            # Read & discard the 16-bit header CRC
            self.fileobj.read(2)
        return dictzipdata

    def write(self, data):
        if self.mode != WRITE:
            import errno
            raise IOError(errno.EBADF, "write() on read-only GzipFile object")
        if self.fileobj is None:
            raise ValueError("write() on closed GzipFile object")
        if len(data) > 0:
            self.offset += len(data)
            if self.dictzip:
                # buffered output
                self.writebuf += data
                self._flush_write_buffer(self.buffersize)
            else:
                # unbuffered output
                self._add_write_data(data)

    def _add_write_data(self, data):
        self.size = self.size + len(data)
        self.crc = zlib.crc32(data, self.crc)
        self.fileobj.write(self.compress.compress(data))

    def _flush_write_buffer(self, buffersize):
        while len(self.writebuf) > buffersize:
            towrite = min(self.buffersize, len(self.writebuf))
            self._write_member(self.writebuf[:towrite])
            self.writebuf = self.writebuf[towrite:]

    def _write_member(self, bufdata):
        # writes complete gzip member (including dictzip table) at once
        self._write_gzip_header(len(bufdata))
        start_of_block = self.fileobj.tell()
        block_sizes = []
        for chunkstart in range(0, len(bufdata), self.chunksize):
            # output the dictzip chunks
            self._add_write_data(
                bufdata[chunkstart: min(chunkstart + self.chunksize, len(bufdata))])
            # flush the compressor and store block size
            self.fileobj.write(self.compress.flush(zlib.Z_FULL_FLUSH))
            current_pos = self.fileobj.tell()
            block_sizes.append(current_pos - start_of_block)
            start_of_block = current_pos
        # finish this member
        self._endmember()
        # seek back and write block sizes
        current_pos = self.fileobj.tell()
        self.fileobj.seek(self.chunktablepos)
        for block_size in block_sizes:
            self.fileobj.write(chr(block_size % 256))
            self.fileobj.write(chr(block_size // 256))
        # return to previous position
        self.fileobj.seek(current_pos)
        # initialize - with no filename - for next member
        self._init_write(None)

    def read(self, size=-1, _block_read_size=None):
        if self.mode != READ:
            import errno
            raise IOError(errno.EBADF, "read() on write-only GzipFile object")

        if self.extrasize <= 0 and self.fileobj is None:
            return ''

        if not _block_read_size:
            # start small, in case the compression factor is high
            _block_read_size = 512
        else:
            _block_read_size = min(_block_read_size, 2048)
        if size < 0:                 # get the whole thing
            try:
                while True:
                    self._read(_block_read_size)
                    _block_read_size = min(
                        _block_read_size * 2, size - self.extrasize + 20)
            except EOFError:
                size = self.extrasize
        else:               # just get some more of it
            try:
                while size > self.extrasize:
                    self._read(_block_read_size)
                    _block_read_size = min(
                        _block_read_size * 2, size - self.extrasize + 20)
            except EOFError:
                if size > self.extrasize:
                    size = self.extrasize

        chunk = self.extrabuf[:size]
        self.extrabuf = self.extrabuf[size:]
        self.extrasize = self.extrasize - size

        self.offset += size
        return chunk

    def _unread(self, buf):
        self.extrabuf = buf + self.extrabuf
        self.extrasize = len(buf) + self.extrasize
        self.offset -= len(buf)

    def _iseof(self):
        pos = self.fileobj.tell()   # Save current position
        self.fileobj.seek(0, 2)     # Seek to end of file
        if pos == self.fileobj.tell():
            raise EOFError("Reached EOF")
        else:
            self.fileobj.seek(pos)  # Return to original position

    def _read(self, size):
        if self.fileobj is None:
            raise EOFError("Reached EOF")

        if self._new_member:
            # If the _new_member flag is set, we have to
            # jump to the next member, if there is one.
            #
            # First, check if we're at the end of the file;
            # if so, it's time to stop; no more members to read.
            self._iseof()

            self._init_read()
            self._read_gzip_header()
            self.decompress = zlib.decompressobj(-zlib.MAX_WBITS)
            self._new_member = False

        # Read a chunk of data from the file
        buf = self.fileobj.read(size)

        # If the EOF has been reached, flush the decompression object
        # and mark this object as finished.

        if buf == "":
            uncompress = self.decompress.flush()
            self._read_eof()
            self._add_read_data(uncompress)
            raise EOFError('Reached EOF')

        uncompress = self.decompress.decompress(buf)
        self._add_read_data(uncompress)

        if self.decompress.unused_data != "":
            # Ending case: we've come to the end of a member in the file,
            # so seek back to the start of the unused data, finish up
            # this member, and read a new gzip header.
            # (The number of bytes to seek back is the length of the unused
            # data, minus 8 because _read_eof() will rewind a further 8 bytes)
            self.fileobj.seek(-len(self.decompress.unused_data) + 8, 1)

            # Check the CRC and file size, and set the flag so we read
            # a new member on the next call
            self._read_eof()
            self._new_member = True

    def _add_read_data(self, data):
        self.crc = zlib.crc32(data, self.crc)
        self.extrabuf = self.extrabuf + data
        self.extrasize = self.extrasize + len(data)
        self.size = self.size + len(data)

    def _read_eof(self):
        # We've read to the end of the file, so we have to rewind in order
        # to reread the 8 bytes containing the CRC and the file size.
        # We check the that the computed CRC and size of the
        # uncompressed data matches the stored values.  Note that the size
        # stored is the true file size mod 2**32.
        if self.dictzip:
            # we don't keep crc values for dictzips
            return
        self.fileobj.seek(-8, 1)
        crc32 = read32(self.fileobj)
        isize = U32(read32(self.fileobj))   # may exceed 2GB
        if U32(crc32) != U32(self.crc):
            raise IOError("CRC check failed")
        elif isize != LOWU32(self.size):
            raise IOError("Incorrect length of data produced")

    def _endmember(self):
        if self.mode == WRITE:
            self.fileobj.write(self.compress.flush())  # unbuffered output
            write32(self.fileobj, self.crc)
            # self.size may exceed 2GB, or even 4GB
            write32u(self.fileobj, LOWU32(self.size))

    def close(self):
        if self.mode == WRITE:
            if self.dictzip:
                # writes buffers, and finishes member
                self._flush_write_buffer(0)
            else:
                self._endmember()
            self.fileobj = None
        elif self.mode == READ:
            self.fileobj = None
        if self.myfileobj:
            self.myfileobj.close()
            self.myfileobj = None

    def __del__(self):
        try:
            if (self.myfileobj is None and
                    self.fileobj is None):
                return
        except AttributeError:
            return
        self.close()

    def flush(self):
        self.fileobj.flush()

    def fileno(self):
        """Invoke the underlying file object's fileno() method.

        This will raise AttributeError if the underlying file object
        doesn't support fileno().
        """
        return self.fileobj.fileno()

    def isatty(self):
        return False

    def tell(self):
        return self.offset

    def rewind(self):
        '''Return the uncompressed stream file position indicator to the
        beginning of the file'''
        if self.mode != READ:
            raise IOError("Can't rewind in write mode")
        self.fileobj.seek(0)
        self._new_member = True
        self.extrabuf = ""
        self.extrasize = 0
        self.offset = 0

    def seek(self, offset, whence=SEEK_SET):
        if whence == SEEK_CUR:
            offset += self.offset
        if self.mode == WRITE:
            if whence == SEEK_END:
                raise IOError('SEEK_END not supported in write mode')
            if offset < self.offset:
                raise IOError('Negative seek in write mode')
            count = offset - self.offset
            for i in range(count // 1024):
                self.write(1024 * '\0')
            self.write((count % 1024) * '\0')
        elif self.mode == READ:
            readahead = None
            if self.dictzip:
                if whence == SEEK_END:
                    offset += self.uncompressed_length
                count = offset - self.offset
                if count >= 0 and count < 32768:
                    # small forward seeks
                    pass
                else:
                    # use table of flush points
                    member = max(
                        0, bisect.bisect_right(self.memberoffset, offset) - 1)
                    # Start of member's uncompressed stream
                    memberoffset = self.memberoffset[member]
                    # Chunk length
                    chlen = self.memberchlen[member]
                    flushpoints = self.memberflushpoints[member]
                    # Exclude the last entry, pointing beyond data stream
                    idx = min(
                        (offset - memberoffset) // chlen, len(flushpoints) - 2)
                    readahead = flushpoints[
                        min(idx + 2, len(flushpoints) - 1)] - flushpoints[idx]
                    # Seek to start of chunk
                    self.fileobj.seek(flushpoints[idx])
                    # Offset of chunk in uncompressed stream
                    self.offset = memberoffset + idx * chlen
                    # Size is relative to member (ignored)
                    self.size = idx * chlen
                    self.extrabuf = ""
                    self.extrasize = 0
                    # CRC is invalid (ignored)
                    self.crc = zlib.crc32("")
                    self.decompress = zlib.decompressobj(-zlib.MAX_WBITS)
                    # Reset possible EOF
                    self._new_member = False
                    count = offset - self.offset
            else:
                if whence == SEEK_END:
                    raise IOError(
                        "SEEK_END only supported on gzip files with a random-access header")
                count = offset - self.offset
                if count < 0:
                    # for backward seek, rewind and do positive seek.
                    # If self.buffersize is set, take this as a hint that a dictzip file
                    # was expected, and refuse to do the (very inefficient)
                    # seek
                    if self.buffersize:
                        raise IOError(
                            "Negative seek on non-dictzip gzip files attempted")
                    self.rewind()
                    count = offset
            # read away unwanted bytes.  A 32K block size appears to be most
            # efficient
            for i in range(count // 32768):
                self.read(32768, _block_read_size=readahead)
            self.read(count % 32768, _block_read_size=readahead)

    def readline(self, size=-1):
        if size < 0:
            size = sys.maxsize
        bufs = []
        readsize = min(100, size)    # Read from the file in small chunks
        while True:
            if size == 0:
                return "".join(bufs)  # Return resulting line

            c = self.read(readsize)
            i = c.find('\n')
            if size is not None:
                # We set i=size to break out of the loop under two
                # conditions: 1) there's no newline, and the chunk is
                # larger than size, or 2) there is a newline, but the
                # resulting line would be longer than 'size'.
                if i == -1 and len(c) > size:
                    i = size - 1
                elif size <= i:
                    i = size - 1

            if i >= 0 or c == '':
                bufs.append(c[:i + 1])    # Add portion of last chunk
                self._unread(c[i + 1:])   # Push back rest of chunk
                return ''.join(bufs)    # Return resulting line

            # Append chunk to list, decrease 'size',
            bufs.append(c)
            size = size - len(c)
            readsize = min(size, readsize * 2)

    def readlines(self, sizehint=0):
        # Negative numbers result in reading all the lines
        if sizehint <= 0:
            sizehint = sys.maxsize
        L = []
        while sizehint > 0:
            line = self.readline()
            if line == "":
                break
            L.append(line)
            sizehint = sizehint - len(line)

        return L

    def writelines(self, L):
        for line in L:
            self.write(line)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if line:
            return line
        else:
            raise StopIteration


def _test():
    # Act like gzip; with -d, act like gunzip; with -D, act like dictzip
    # The input file is not deleted, however, nor are any other gzip
    # options or features supported.
    args = sys.argv[1:]
    decompress = args and args[0] == "-d"
    dictzip = args and args[0] == "-D"
    if decompress or dictzip:
        args = args[1:]
    if not args:
        args = ["-"]
    for arg in args:
        if decompress:
            if arg == "-":
                f = GzipFile(filename="", mode="rb", fileobj=sys.stdin)
                g = sys.stdout
            else:
                if arg[-3:] != ".gz" and arg[-3:] != ".dz":
                    print("filename doesn't end in .gz or .dz:", repr(arg))
                    continue
                f = open(arg, "rb")
                g = builtins.open(arg[:-3], "wb")
        else:
            if dictzip:
                buffersize = 1000000
                ext = ".dz"
            else:
                buffersize = None
                ext = ".gz"
            if arg == "-":
                f = sys.stdin
                g = GzipFile(
                    filename="", mode="wb", fileobj=sys.stdout, chunksize=1000, buffersize=buffersize)
            else:
                f = builtins.open(arg, "rb")
                g = open(
                    arg + ext, "wb", chunksize=1000, buffersize=buffersize)
        blocksize = 1024
        if False:
            while True:
                chunk = f.read(blocksize)
                if not chunk:
                    break
                g.write(chunk)
        else:
            # test the random access code
            ptr = 0
            while True:
                f.seek(0)
                f.seek(ptr)
                chunk = f.read(blocksize)
                if not chunk:
                    break
                g.write(chunk)
                ptr += blocksize
        if g is not sys.stdout:
            g.close()
        if f is not sys.stdin:
            f.close()

if __name__ == '__main__':
    _test()
