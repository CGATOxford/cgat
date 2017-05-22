'''
Mali.py - Tools for multiple alignments 
=======================================

:Tags: Python

Code
----

'''

import re
import string
import math
import copy
import random


class AlignedString:

    mGapChars = ("-", ".")
    mGapChar = "-"

    def __init__(self, identifier, fr, to, s, is_mangled=False):
        self.mId = identifier
        self.mFrom = fr
        self.mTo = to
        self.mString = s
        self.mIsMangled = is_mangled

    def __len__(self):
        return len(self.mString)

    def maskColumn(self, column, mask_char="x"):
        self.maskColumns(self, [column, ], mask_char=mask_char)

    def maskColumns(self, columns, mask_char="x"):
        s = list(self.mString)
        for c in columns:
            if s[c] not in self.mGapChars:
                s[c] = mask_char
        self.mString = "".join(s)

    def mapColumns(self, columns, map_function):
        s = list(self.mString)
        for c in columns:
            if s[c] not in self.mGapChars:
                s[c] = map_function(s[c])
        self.mString = "".join(s)

    def takeColumns(self, columns):
        """take columns.

        not implemented yet: updating of residue numbers.
        """
        s = []

        for c in columns:
            s.append(self.mString[c])

        self.mString = "".join(s)

    def truncate(self, first, last):
        """truncate aligned string.

        not implemented yet: updating of residue numbers.
        """
        self.mString = self.mString[first:last]

    def insertColumns(self, position, num_columns, char="-"):
        """insert num_columns columns at position."""
        self.mString = self.mString[:position] + \
            char * num_columns + self.mString[position:]

    def getIdentifier(self):
        '''return demangled identifier.'''
        if self.mIsMangled:
            return self.mId.rsplit("_", 1)[0]
        else:
            return self.mId

    def getResidueNumber(self, position):
        x = self.mFrom
        for c in self.mString[:position]:
            if c not in self.mGapChars:
                x += 1

        return x

    def getSequence(self):
        """return sequence without gaps."""
        return re.sub("[%s]" % "".join(self.mGapChars), "", self.mString)

    def getSegments(self, transitions, gap_chars="-."):
        """return segments of alignment according to
        transitions.

        Transitions is a list of residue numbers.
        """
        segments = []
        pos = self.mFrom

        first_x = 0
        for x in range(len(self.mString)):

            char = self.mString[x]

            if char not in gap_chars:
                if pos in transitions:
                    if first_x != x:
                        segments.append((first_x, x))
                    first_x = x
                pos += 1

        # add last segment, unless it has already been added
        # because it was part of the transitions
        if first_x != x:
            segments.append((first_x, len(self.mString)))
        elif pos in transitions:
            segments.append((first_x, len(self.mString)))

        return segments

    def threadSequence(self, new_sequence, map_old2new):
        """thread a new sequence onto this sequence replacing
        characters in this sequence with those in new sequence.

        ``map_old2new`` lists the positions that should be replaced.
        It maps non-gap positions in old to those in new. Non-gap
        positions than are not mapped are set to :attr:`gapchar`.
        """

        s = list(self.mString)
        c = 0
        for x in range(len(s)):
            m = map_old2new.mapRowToCol(c)
            if m >= 0:
                s[c] = new_sequence[m]
            else:
                s[c] = self.mGapChar
            if s[c] not in self.mGapChars:
                c += 1

        self.mString = "".join(s)


class MaliData:

    def __init__(self, line, gap_chars="-.", mask_chars="Nn"):

        self.mNGaps = len(re.sub("[^%s]" % gap_chars, "", line))
        self.mNMasked = len(re.sub("[^%s]" % mask_chars, "", line))
        self.mNAll = len(line)
        self.mNChars = self.mNAll - self.mNMasked - self.mNGaps

    def __str__(self):
        return "%i\t%i\t%i\t%i" % (self.mNAll, self.mNChars,
                                   self.mNGaps, self.mNMasked)


class Mali:

    mGapChars = ("-", ".")
    mGapPattern = re.compile("[.-]")
    # character used for gaps
    mGapChar = "-"
    # admissable gap characters
    mGapChars = ("-", ".")
    mMaskChar = "X"

    def __init__(self):

        self.mIdentifiers = []
        self.mMali = {}
        self.mLength = 0
        self.mAnnotations = {}
        self.mName = None

        # set to false, if ranges shall not be output.
        self.mWriteRanges = True

    def __contains__(self, key):
        return key in self.mMali

    def __getitem__(self, key):
        return self.mMali[key].mString

    def __len__(self):
        return len(self.mIdentifiers)

    def __delitem__(self, key):
        del self.mMali[key]
        self.mIdentifiers.remove(key)

    def getClone(self):
        return copy.deepcopy(self)

    def iteritems(self):
        return iter(self.mMali.items())

    def items(self):
        return list(self.mMali.items())

    def values(self):
        return list(self.mMali.values())

    def keys(self):
        return list(self.mMali.keys())

    def getName(self):
        return self.mName

    def setName(self, name):
        self.mName = name

    def getIdentifiers(self):
        return self.mIdentifiers

    def getLength(self):
        """deprecated."""
        return self.getNumSequences()

    def getNumSequences(self):
        return self.__len__()

    def getNumColumns(self):
        if self.mIdentifiers:
            return len(self.mMali[self.mIdentifiers[0]])
        else:
            return 0

    def getWidth(self):
        """deprecated."""
        return self.getNumColumns()

    def rename(self, old_name, new_name):
        """rename an entry."""
        if old_name not in self.mMali:
            raise KeyError("%s not in mali" % old_name)
        self.mIdentifiers[self.mIdentifiers.index(old_name)] = new_name
        self.mMali[new_name] = self.mMali[old_name]
        del self.mMali[old_name]

    def isEmpty(self):
        return len(self.mMali) == 0

    def getSequence(self, key):
        return self.mMali[key]

    def getResidueNumber(self, key, position):
        """return residue number in sequence key at position position."""
        return self.mMali[key].getResidueNumber(position)

    def setSequence(self, key, sequence):
        # TODO: updating of residue numbers
        self.mMali[key].mString = sequence

    def getEntry(self, key):
        return self.mMali[key]

    def addEntry(self, s):
        """add an aligned string object."""
        if s.mId in list(self.mMali.keys()):
            raise KeyError("id %s already in mali" % s.mId)

        self.mIdentifiers.append(s.mId)
        self.mMali[s.mId] = s

    def addSequence(self, id, fr, to, sequence):

        if to < 0:
            to = self.countCharacters(sequence)
        s = AlignedString(id, fr, to, sequence)
        self.addEntry(s)

    def countCharacters(self, row):
        return len(row) - len(self.mGapPattern.findall(row))

    def deleteEntry(self, identifier):
        if identifier not in self.mMali:
            raise KeyError("identifier %s not in mali." % identifier)
        del self.mMali[identifier]
        self.mIdentifiers.remove(identifier)

    def getColumns(self):
        """return mali in column orientation."""
        args = [self.mMali[x].mString for x in self.mIdentifiers]
        return ["".join(x) for x in zip(*args)]

    def getConsensus(self, mark_with_gaps=False):
        """return consensus string.

        The consensus string returns the most frequent character per column
        that is not a gap.
        If mark_with_gaps is set to True, positions with any gap characater are
        set to gaps.
        """
        columns = self.getColumns()
        seq = []
        for x in range(len(columns)):
            s = columns[x]
            counts = [(a, s.count(a))
                      for a in set(list(s)).difference(set(self.mGapChar))]
            if mark_with_gaps and self.mGapChar in columns[x]:
                seq.append(self.mGapChar)
            else:
                counts.sort(key=lambda x: -x[1])
                seq.append(counts[0][0])
        return "".join(seq)

    def readFromFile(self, infile, format="fasta"):
        """read multiple alignment from file in various format."""

        self.mMali = {}
        self.mIdentifiers = []

        pattern_parse_ranges = re.compile("(\S+)/(\d+)-(\d+)")

        # read profiles - a profile possibly consists of several
        # entries per file so treat it differently
        if format.lower() == "profile":

            while 1:
                line = infile.readline()
                if not line:
                    return False
                if line[0] != "#":
                    break

            if line[0] != ">":
                raise "expected '>' at as first character in line %s" % line

            try:
                self.mName, length, width = re.match(
                    ">profile=(\S+) length=(\d+) width=(\d+)", line).groups()
            except AttributeError:
                raise "could not parse header line %s" % line
            width = int(width)
            for x in range(0, width):
                id = "seq%i" % x
                self.mIdentifiers.append(id)
                line = infile.readline()
                if not line:
                    raise "expected %i sequences, only got %i" % (width, x)
                self.mMali[id] = AlignedString(
                    id, 0, self.countCharacters(line[:-1]), line[:-1])
            return True

        if isinstance(infile, list) or \
           isinstance(infile, tuple):
            lines = infile
        else:
            lines = infile.readlines()

        if format not in ("stockholm"):
            # save comments
            self.mComments = [x for x in lines if x[0] == "#"]
            lines = [x for x in lines if x[0] != "#"]
        else:
            self.mComments = []

        # remove empty lines
        lines = [x for x in lines if x.strip()]
        if not lines:
            raise AttributeError("empty alignment")

        def getId(id, s):
            x = pattern_parse_ranges.match(id)
            if x:
                id, fr, to = x.groups()
                fr, to = int(fr) - 1, int(to)
            else:
                fr, to = 0, self.countCharacters(s)
                self.mWriteRanges = False

            return id, fr, to

        #######################################################################
        if format.lower() == "plain":

            for line in lines:
                if not line.strip():
                    continue
                data = line[:-1].split("\t")
                id = data[3]
                xid = id
                x = 0
                while xid in self.mMali:
                    xid = id + "-" + str(x)
                    x += 1

                self.addEntry(
                    AlignedString(xid,
                                  int(data[0]) - 1,
                                  int(data[2]),
                                  data[1]))

        #######################################################################
        elif format.lower() == "fasta":
            pattern_identifier = "\S+"
            id = None
            fragments = []
            for line in lines:
                if line[0] == ">":
                    if id:
                        s = re.sub("\s", "", string.join(fragments, ""))
                        id, fr, to = getId(id, s)
                        self.addEntry(AlignedString(id, fr, to, s))

                    id = re.search(
                        "^(%s)" % pattern_identifier, line[1:-1]).group(0)
                    fragments = []
                    continue
                fragments.append(line[:-1])

            s = re.sub("\s", "", string.join(fragments, ""))
            id, fr, to = getId(id, s)
            self.addEntry(AlignedString(id, fr, to, s))

        #######################################################################
        elif format.lower() == "phylip":
            nsequences, nchars = re.split("\s+", lines[0][:-1].strip())
            nsequences = int(nsequences)
            for line in lines[1:]:
                l = line[:-1].strip()
                if not l:
                    continue
                id, sequence = re.match("(\S+)\s+(.*)", l).groups()
                sequence = re.sub("\s", "", sequence)
                if id not in self.mMali:
                    self.mIdentifiers.append(id)
                    self.mMali[id] = []

                self.mMali[id].append(sequence)

            for id, frags in list(self.mMali.items()):
                s = "".join(frags)
                fr, to = 0, self.countCharacters(s)
                self.mMali[id] = AlignedString(id, fr, to, s)

        #######################################################################
        elif format.lower() == "clustal":
            # skip header line
            del lines[0]
            fragments = {}

            # prune lines
            lines = [x.strip() for x in lines]
            # remove empty lines
            lines = [x for x in lines if len(x[:-1]) > 0]

            for line in lines:
                # remove consensus lines
                if line[0] in ("*", ":"):
                    continue

                data = re.split("\s+", line)
                if len(data) != 2:
                    raise ValueError("parsing error in line %s" % line)

                id, fragment = data
                if id not in fragments:
                    fragments[id] = []
                    self.mIdentifiers.append(id)

                fragments[id].append(fragment)

            for id, f in list(fragments.items()):
                s = re.sub("\s", "", string.join(f, ""))
                self.mMali[id] = AlignedString(
                    id, 0, self.countCharacters(s), s)

        elif format.lower() == "stockholm":
            # skip header line
            assert lines[0].startswith(
                "# STOCKHOLM"), "file is not in stockholm format"
            del lines[0]
            fragments = {}
            annotations = {}
            # prune lines
            lines = [x.strip() for x in lines]
            # remove empty lines
            lines = [x for x in lines if len(x[:-1]) > 0]

            for line in lines:
                data = re.split("\s+", line)

                if data[0] == "//":
                    break

                if line[0] == '#':
                    if data[0] == "#=GC":
                        id, fragment = data[1:3]
                    else:
                        self.mComments.append(line)
                        continue
                    if id not in annotations:
                        annotations[id] = []
                    annotations[id].append(fragment)
                else:

                    if len(data) > 2:
                        raise ValueError("parsing error in line %s" % line)
                    elif len(data) == 1:
                        # treat empty alignments/lines
                        id = data[0]
                        fragment = ""
                    else:
                        id, fragment = data

                    if id not in fragments:
                        fragments[id] = []
                        self.mIdentifiers.append(id)

                    fragments[id].append(fragment)

            n = []
            for id in self.mIdentifiers:
                f = fragments[id]

                s = re.sub("\s", "", string.join(f, ""))
                x = pattern_parse_ranges.match(id)
                if x:
                    id, fr, to = x.groups()
                    fr, to = int(fr) - 1, int(to)
                else:
                    fr, to = 0, self.countCharacters(s)

                n.append(id)
                self.mMali[id] = AlignedString(id, fr, to, s)
            self.mIdentifiers = n

            for id, f in list(annotations.items()):
                s = re.sub("\s", "", string.join(f, ""))
                annotations[id] = s
            self.mAnnotations = annotations
        else:
            raise "unknown alignment format %s" % format

        if len(self.mMali) == 0:
            self.mLength = 0
        else:
            self.mLength = min(
                [len(x.mString) for x in list(self.mMali.values())])

    def writeToFile(self, outfile, write_ranges=True, format="plain",
                    options=None):
        """write alignment to file.

        If options is given, these lines are output into the multiple
        alignment.

        """
        if format == "plain-fasta":
            format = "fasta"
            write_ranges = False

        write_ranges = write_ranges and self.mWriteRanges

        if format == "plain":
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                outfile.write("%i\t%s\t%i\t%s\n" % (
                    m.mFrom + 1, m.mString, m.mTo, m.getIdentifier()))

        elif format == "fasta":
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                if write_ranges:
                    outfile.write(
                        ">%s/%i-%i\n%s\n" % (m.getIdentifier(), m.mFrom + 1, m.mTo, m.mString))
                else:
                    outfile.write(">%s\n%s\n" % (identifier, m.mString))

        elif format == "stockholm":
            outfile.write("# STOCKHOLM 1.0\n")

            if options:
                for o in options:
                    outfile.write("%s\n" % o)

            # calculate offset:
            max_l = 0
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                id = m.getIdentifier()
                # tab does not work as separator
                if m.mTo and write_ranges:
                    x = "%s/%i-%i" % (id, m.mFrom + 1, m.mTo)
                else:
                    x = "%s" % (id)
                max_l = max(max_l, len(x))
            for identifier in list(self.mAnnotations.keys()):
                x = "#=GC %s" % identifier
                max_l = max(max_l, len(x))

            format = "%-" + str(max_l) + "s  %s\n"
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                id = m.getIdentifier()
                # tab does not work as separator
                if m.mTo and write_ranges:
                    x = "%s/%i-%i" % (id, m.mFrom + 1, m.mTo)
                else:
                    x = "%s" % (id)

                outfile.write(format % (x, m.mString))

            for identifier, value in list(self.mAnnotations.items()):
                x = "#=GC %s" % identifier
                outfile.write(format % (x, value))

            outfile.write("//\n")

        elif format == "phylip":

            outfile.write("%i %i\n" % (self.getLength(), self.getWidth()))

            for identifier in self.mIdentifiers:
                outfile.write("%s   %s\n" %
                              (identifier, self.mMali[identifier].mString))

        elif format.lower() == "profile":
            if self.mName:
                name = self.mName
            else:
                name = ",".join(self.mIdentifiers)

            outfile.write(">profile=%s length=%i width=%i\n" %
                          (name, self.getWidth(), self.getLength()))

            for identifier in self.mIdentifiers:
                outfile.write("%s\n" % (self.mMali[identifier].mString))

        elif format == "nexus":
            # nexus formatted output - MrBayes conformant.
            outfile.write("#NEXUS\n")
            outfile.write("begin data;\n")
            outfile.write("  dimensions ntax=%i nchar=%i;\n" %
                          (self.getLength(), self.getWidth()))
            outfile.write(
                "  format datatype=dna interleave=no gap=%s;\n" % (self.mGapChar))
            outfile.write("  matrix\n")
            max_len = max([len(x) for x in self.mIdentifiers])
            format = "  %-" + str(max_len) + "s %s\n"
            for identifier in self.mIdentifiers:
                outfile.write(
                    format % (identifier, self.mMali[identifier].mString))
            outfile.write("  ;\n")
            outfile.write("end;\n")

        else:
            raise "unknown alignment format %s" % format

    def removeUnalignedEnds(self):
        """remove unaligned ends in the multiple alignment.

        unaligned ends correspond to lower-case characters.
        """
        pattern_start = re.compile("^([- .a-z]+)")
        pattern_unaligned = re.compile("[a-z]")

        for s in list(self.mMali.values()):

            first = pattern_start.match(s.mString)
            if first:
                first = first.groups()[0]
                nchars = len(pattern_unaligned.findall(first))
                s.mFrom += nchars
                s.mString = self.mGapChar * len(first) + s.mString[len(first):]

            # search from the back end by reversing. This is much faster than
            # using $ from the back.
            last = pattern_start.match(s.mString[::-1])
            if last:
                last = last.groups()[0]
                nchars = len(pattern_unaligned.findall(last))
                s.mTo -= nchars
                l = len(s) - len(last)
                s.mString = s.mString[:l] + self.mGapChar * l

    def upperCase(self):
        """set all characters to upper case."""
        for k, s in list(self.mMali.items()):
            s.mString = s.mString.upper()

    def lowerCase(self):
        """set all characters to lower case."""
        for k, s in list(self.mMali.items()):
            s.mString = s.mString.lower()

    def removeEndGaps(self):
        """remove end gaps.

        end gaps do not include any characters and thus
        the alignment coordinates won't change.
        """

        pattern_start_gaps = re.compile("^([- ]+)")

        min_from = self.mLength
        max_to = 0

        for s in list(self.mMali.values()):

            first = pattern_start_gaps.match(s.mString)
            if first:
                first = first.groups()[0]
                min_from = min(min_from, len(first))

            # search from the back end by reversing. This is much faster than
            # using $ from the back.
            last = pattern_start_gaps.search(s.mString[::-1])
            if last:
                last = last.groups()[0]
                max_to = max(max_to, len(s) - len(last))

        for s in list(self.mMali.values()):
            s.mString = s.mString[min_from:max_to]

        self.mLength = min([x.mString for x in list(self.mMali.values())])

    def insertColumns(self, position, num_gaps, keep_fixed=None, char="-"):
        """insert gaps at position into multiple alignment.

        if keep_constant is a list of identifiers, those are kept constant,
        instead, gaps are added to the end.
        """

        last_pos = min(self.getWidth(), position + num_gaps)
        for id, seq in list(self.mMali.items()):
            if keep_fixed and id in keep_fixed:
                seq.insertColumns(last_pos, num_gaps, char)
            else:
                seq.insertColumns(position, num_gaps, char)

    def removeGaps(self,
                   allowed_gaps=0,
                   minimum_gaps=1,
                   frame=1):
        """remove gappy columns.

        allowed_gaps: number of gaps allowed for column to be kept
        minimum_gaps: number of gaps for column to be removed

        set minimum_gaps to the number of sequences to remove columns
        with all gaps.

        If frame is > 1 (3 most likely), then a whole codon is removed
        as soon as there is one column to be removed.
        """

        self.removePattern(
            match_function=lambda x: x in self.mGapChars,
            allowed_matches=allowed_gaps,
            minimum_matches=minimum_gaps,
            delete_frame=frame)

    def removePattern(self,
                      match_function,
                      allowed_matches=0,
                      minimum_matches=1,
                      delete_frame=1,
                      search_frame=1):
        """remove columns (or group of columns), that match a certain pattern.

        allowed_matches: number of matches allowed so that column is still kept
        minimum_matches: number of matches required for column to be removed

        set minimum_matches to the number of sequences to remove columns
        with all gaps.

        Patterns are matches in search_frame. For example, if frame is
        3, whole codons are supplied to match_function.

        delete_frame specifies the frame for deletion. If it is set to 3,
        codons are removed if already one column matches.

        Example: remove all columns that contain at least one stop-codon:

        removePattern( lambda x: x.upper() in ("TAG", "TAA", "TGA"),
                allowed_matches = 0,
                minimum_matches = 1,
                search_frame = 3,
                delete_frame = 3)

        """
        nmatches = [0] * self.getWidth()

        for s in [x.mString for x in list(self.mMali.values())]:
            for x in range(0, len(s), search_frame):
                segment = s[x:x + search_frame]
                if match_function(segment):
                    nmatches[x] += 1

        columns = []
        delete_columns = []

        for x in range(len(nmatches)):
            if nmatches[x] >= allowed_matches and nmatches[x] < minimum_matches:
                columns.append(x)
            else:
                delete_columns.append(x)

        if delete_frame != 1:
            s = set(columns)
            for x in delete_columns:
                start = int(math.floor(float(x) / delete_frame) * delete_frame)
                end = start + delete_frame
                for c in range(start, end):
                    if c in s:
                        s.remove(c)

            columns = list(s)
            columns.sort()

        self.takeColumns(columns)

    def removeEmptySequences(self):
        """remove sequences that are completely empty.
        """
        new_ids = []
        for id in self.mIdentifiers:
            if self.countCharacters(self.mMali[id].mString) == 0:
                del self.mMali[id]
                continue
            new_ids.append(id)
        self.mIdentifiers = new_ids

    def upper(self):
        """convert all characters in mali to uppercase."""

        for s in list(self.mMali.values()):
            s.mString = s.mString.upper()

    def lower(self):
        """convert all characters in mali to lowercase."""

        for s in list(self.mMali.values()):
            s.mString = s.mString.lower()

    def shiftAlignment(self, map_id2offset):
        """shift alignment by offset."""

        for identifier, m in list(self.mMali.items()):
            if identifier in map_id2offset:
                o = map_id2offset[identifier]
                m.mFrom += o
                m.mTo += o

    def markCodons(self, mode="case"):
        """mark codons.
        """
        for identifier, m in list(self.mMali.items()):
            s = m.mString
            if len(s) % 3 != 0:
                raise "sequence %s not divisible by 3" % (m.mId)

            is_upper = True
            sequence = []
            for x in range(0, len(s), 3):
                if is_upper:
                    sequence.append(s[x:x + 3].upper())
                    is_upper = False
                else:
                    sequence.append(s[x:x + 3].lower())
                    is_upper = True

            m.mString = "".join(sequence)

    def markTransitions(self, map_id2transitions, mode="case"):
        """mark transitions in the multiple alignment.

        if mode == case, then upper/lower case is used for the transitions

        Otherwise, a character given by mode is inserted.
        """

        if mode in ("case", "keep-odd", "keep-even"):

            # check, if the whole alignment needs to be masked/marked:
            if "mali" in map_id2transitions:
                transitions = map_id2transitions["mali"]
                for identifier, s in list(self.mMali.items()):
                    new_chars = []
                    is_upper = True
                    is_first = False

                    for c in range(len(s)):

                        if c in transitions:
                            is_first = True
                            if is_upper:
                                is_upper = False
                            else:
                                is_upper = True

                        x = s.mString[c]

                        if mode == "case":
                            if x in string.lowercase:
                                x = self.mMaskChar
                            if is_upper:
                                x = string.upper(x)
                            else:
                                x = string.lower(x)
                        elif mode == "keep-even":
                            if is_upper:
                                x = self.mGapChar
                        elif mode == "keep-odd":
                            if not is_upper:
                                x = self.mGapChar

                        new_chars.append(x)

                    s.mString = "".join(new_chars)

            # now do individual sequences
            for identifier, s in list(self.mMali.items()):
                if identifier not in map_id2transitions:
                    continue

                new_chars = []
                c = s.mFrom

                is_upper = True
                transitions = map_id2transitions[identifier]
                for x in s.mString:
                    if x in self.mGapChars:
                        pass
                    else:
                        if c in map_id2transitions[identifier]:
                            if is_upper:
                                is_upper = False
                            else:
                                is_upper = True
                        c += 1

                        if x in string.lowercase:
                            x = self.mMaskChar

                        if mode == "case":
                            if is_upper:
                                x = string.upper(x)
                            else:
                                x = string.lower(x)
                        elif mode == "keep-even":
                            if is_upper:
                                x = self.mGapChar
                        elif mode == "keep-odd":
                            if not is_upper:
                                x = self.mGapChar

                    new_chars.append(x)

                s.mString = "".join(new_chars)
        else:
            raise ValueError("character insertion not implemented yet")

    def buildColumnMap(self, other, join_field=None):
        """build map of columns in other to this."""

        if not join_field:
            join_field = other.mIdentifiers[0]

        if join_field not in other.mMali or \
           join_field not in self.mMali:
            raise "line %s not in both alignments." % (join_field)

        this_seq = self.mMali[join_field]
        other_seq = other.mMali[join_field]

        if this_seq.mFrom != other_seq.mFrom or \
           this_seq.mTo != other_seq.mTo:
            raise ValueError("residue ranges for sequence %s doe not correspond." % (
                join_field))

        map_this2other = []

        this_seq = this_seq.mString.upper()
        other_seq = other_seq.mString.upper()

        # position in other
        o = 0
        for c in this_seq:
            if c in self.mGapChars:
                map_this2other.append(None)
            else:
                while other_seq[o] != c:
                    o += 1
                map_this2other.append(o)
                o += 1
        return map_this2other

    def shuffle(self, frame=1):
        """shuffle multiple alignment.

        The frame determines the block size for shuffling. Use 3 for codons in
        a multiple alignment without frame-shifts.
        """
        columns = list(range(self.getNumColumns() // frame))
        random.shuffle(columns)
        if frame > 1:
            cc = []
            for x in columns:
                cc += [x * frame + y for y in range(frame)]
            columns = cc
        self.takeColumns(columns)

    def propagateMasks(self, min_chars=1, mask_char="x"):
        """propagate masked characters to all rows of a multiple alignment
        within a column.

        If there is at least min_chars in a mali column, that are masks,
        propagate the masks to all other rows.
        """

        masks_per_column = {}
        for identifier, s in list(self.mMali.items()):
            r = s.mString.lower()
            for x in range(len(r)):
                if r[x] == mask_char:
                    if x not in masks_per_column:
                        masks_per_column[x] = 0
                    masks_per_column[x] += 1

        columns_to_mask = []
        for c, n in list(masks_per_column.items()):
            if n >= min_chars:
                columns_to_mask.append(c)

        columns_to_mask.sort()

        self.maskColumns(columns_to_mask, mask_char=mask_char)

    def propagateTransitions(self, min_chars=1):
        """propagate lower case in a column to all residues.
        """
        columns_to_change = set()
        for identifier, s in list(self.mMali.items()):
            r = s.mString
            for x in range(len(r)):
                if r[x] in string.lowercase:
                    columns_to_change.add(x)

        columns_to_change = list(columns_to_change)
        columns_to_change.sort()

        self.mapColumns(columns_to_change, string.lower)

    def takeColumns(self, columns):
        """restrict alignments to certain columns."""
        for identifier, s in list(self.mMali.items()):
            s.takeColumns(columns)

        for key, anno in list(self.mAnnotations.items()):
            self.mAnnotations[key] = "".join([anno[c] for c in columns])

    def maskColumns(self, columns, mask_char="x"):
        """mask columns in a multiple alignment."""

        for identifier, s in list(self.mMali.items()):
            s.maskColumns(columns, mask_char=mask_char)

    def mapColumns(self, columns, map_function):
        """apply map_function to all residues in columns."""

        for identifier, s in list(self.mMali.items()):
            s.mapColumns(columns, map_function)

    def recount(self, reset_first=False):
        """recount residue in alignments."""
        for id, seq in list(self.mMali.items()):
            if reset_first:
                seq.mFrom = 0
            seq.mTo = seq.mFrom + self.countCharacters(seq.mString)

    def maskColumn(self, column, mask_char="x"):
        """mask a column."""

        for identifier, s in list(self.mMali.items()):
            s.maskColumn(column, mask_char)

    def copyAnnotations(self, other):
        """copy annotations from annother mali."""

        map_this2other = self.buildColumnMap(other)
        ncols = self.getWidth()

        for key, annotation in list(other.mAnnotations.items()):
            a = []
            for x in range(ncols):
                m = map_this2other[x]
                if m is not None:
                    a.append(annotation[m])
                else:
                    a.append(self.mGapChar)

            self.mAnnotations[key] = "".join(a)

    def getAnnotation(self, key):
        """return annotation associated with key."""
        return self.mAnnotations[key]

    def setAnnotation(self, key, value):
        """set annotation associated with key to value."""
        self.mAnnotations[key] = value

    def addAnnotation(self, key, annotation):
        """add annotation."""
        xkey = key
        n = 0
        while xkey in self.mAnnotations:
            xkey = "%s_%i" % (key, n)

        self.mAnnotations[xkey] = re.sub("\s", "", annotation)

    def truncate(self, first, last):
        """truncate alignment within range."""
        for key, value in list(self.mMali.items()):
            value.truncate(first, last)

    def mapIdentifiers(self, map_old2new=None, pattern_identifier="ID%06i"):
        """map identifiers in multiple aligment.

        if map_old2new is not given, a new map is created (map_new2old)
        """
        new_identifiers = []
        if map_old2new is not None:
            for id in self.mIdentifiers:
                new_id = map_old2new[id]
                new_identifiers.append(new_id)

                entry = self.mMali[id]
                del self.mMali[id]
                # entry.mId = new_id
                self.mMali[new_id] = entry
        else:
            map_old2new = {}
            n = 1
            for id in self.mIdentifiers:
                new_id = pattern_identifier % n
                new_identifiers.append(new_id)

                entry = self.mMali[id]
                del self.mMali[id]
                # entry.mId = new_id
                self.mMali[new_id] = entry

                map_old2new[new_id] = id
                n += 1

        self.mIdentifiers = new_identifiers

        return map_old2new

    def getAlphabet(self):
        """get alphabet from the multiple alignment.

        Alphabet is "na", if more than 90% of characaters are "actgxn",
        otherwise it is "aa".
        """
        s = "".join([x.mString for x in list(self.values())]).lower()
        s = re.sub("[%s]" % "".join(self.mGapChars), "", s)
        ss = re.sub("[acgtxn]", "", s)
        if float(len(ss)) < (len(s) * 0.1):
            return "na"
        else:
            return "aa"

    def checkLength(self):
        """check lengths of aligned strings.

        Return false if they are inconsistent."""

        l = None
        for s in list(self.mMali.values()):
            if l is None:
                l = len(s.mString)
            else:
                if l != len(s.mString):
                    return False
        return True

    def clipByAnnotation(self, key, chars=""):
        """restrict alignment to positions where annotation identified by key
        in chars.

        if chars is empty, nothing is clipped.

        """
        if key not in self.mAnnotations:
            return
        if chars == "":
            return
        if self.getNumColumns() == 0:
            return

        anno = self.mAnnotations[key]
        columns = []

        for x in range(len(anno)):
            if anno[x] in chars:
                columns.append(x)

        self.takeColumns(columns)

    def apply(self, f):
        """apply function f to every row in the multiple alignment."""
        for s in list(self.mMali.values()):
            f(s)

    def filter(self, f):
        '''filter multiple alignment using function *f*.

        The function *f* should return True for entries that should be
        kept and False for those that should be removed.

        '''
        ids, vals = [], []
        for id in self.mIdentifiers:
            s = self.mMali[id]
            if f(s):
                vals.append(s)
                ids.append(id)
        self.mIdentifiers = ids
        self.mMali = vals


class SequenceCollection(Mali):

    """reads in a sequence collection, but permits several entries per id.

    Note that this might cause problems with interleaved formats like
    phylips or clustal.

    This mapping is achieved by appending a numeric suffix. The suffix
    is retained for the life-time of the object, but not output to a
    file.

    """

    def __init__(self):
        Mali.__init__(self)

    def addEntry(self, s):
        """add an aligned string object."""
        id = s.mId
        if id in list(self.mMali.keys()):
            x = 1
            while "%s_%i" % (id, x) in list(self.mMali.keys()):
                x += 1

            id = "%s_%i" % (id, x)
        self.mIdentifiers.append(id)
        s.mIsMangled = True

        self.mMali[id] = s

    def readFromFile(self, infile, format="fasta"):
        """read multiple alignment from file in various format."""

        self.mMali = {}
        self.mIdentifiers = []

        pattern_parse_ranges = re.compile("(\S+)/(\d+)-(\d+)")

        if isinstance(infile, list) or \
           isinstance(infile, tuple):
            lines = infile
        else:
            lines = infile.readlines()

        if format not in ("stockholm"):
            # save comments
            self.mComments = [x for x in lines if x[0] == "#"]
            lines = [x for x in lines if x[0] != "#"]
        else:
            self.mComments = []

        def getId(id, s):
            x = pattern_parse_ranges.match(id)
            if x:
                id, fr, to = x.groups()
                fr, to = int(fr) - 1, int(to)
            else:
                fr, to = 0, self.countCharacters(s)
                self.mWriteRanges = False

            return id, fr, to

        #######################################################################
        if format.lower() == "plain":

            for line in lines:
                if not line.strip():
                    continue
                data = line[:-1].split("\t")
                id = data[3]
                xid = id
                x = 0
                while xid in self.mMali:
                    xid = id + "-" + str(x)
                    x += 1

                self.addEntry(
                    AlignedString(xid,
                                  int(data[0]) - 1,
                                  int(data[2]),
                                  data[1]))

        #######################################################################
        elif format.lower() == "fasta":
            pattern_identifier = "\S+"
            id = None
            fragments = []
            for line in lines:
                if line[0] == ">":
                    if id:
                        s = re.sub("\s", "", string.join(fragments, ""))
                        id, fr, to = getId(id, s)
                        self.addEntry(AlignedString(id, fr, to, s))

                    id = re.search(
                        "^(%s)" % pattern_identifier, line[1:-1]).group(0)
                    fragments = []
                    continue
                fragments.append(line[:-1])

            s = re.sub("\s", "", string.join(fragments, ""))
            id, fr, to = getId(id, s)
            self.addEntry(AlignedString(id, fr, to, s))

        #######################################################################
        elif format.lower() == "phylip":
            raise ValueError("phylip not implemented")

        #######################################################################
        elif format.lower() == "clustal":
            raise ValueError("clustal not implemented")

        elif format.lower() == "stockholm":
            raise ValueError("stockholm not implemented")

        if len(self.mMali) == 0:
            self.mLength = 0
        else:
            self.mLength = min(
                [len(x.mString) for x in list(self.mMali.values())])


def convertMali2Alignlib(mali):
    '''convert a multiple alignment of type :class:`Mali`
    into an alignlib_lite.py_multiple alignment object.
    '''

    import alignlib_lite
    m = alignlib_lite.py_makeMultipleAlignment()
    for identifier in mali.getIdentifiers():
        a = alignlib_lite.py_makeAlignatum(mali[identifier])
        m.add(a)
    return m


def convertAlignlib2Mali(mali, identifiers=None, seqs=None):
    """convert a multiple alignment into an alignlib_lite.py_multiple
    alignment object."""
    m = Mali()

    if not identifiers:
        identifiers = ["%i" % x for x in range(mali.getNumSequences())]

    if seqs is None:
        # old style MultipleAlignment
        for x in range(mali.getNumSequences()):
            a = mali.getRow(x)
            m.addSequence(
                identifiers[x], a.getFrom(), a.getTo(), a.getString())
    else:
        import alignlib_lite
        output = alignlib_lite.py_MultAlignmentFormatPlain(mali, seqs)
        for x in range(mali.getNumSequences()):
            a = output.mData[x]
            m.addSequence(
                identifiers[x], a.getFrom(), a.getTo(), a.getString())

    return m
