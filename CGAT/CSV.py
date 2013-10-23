"""
CSV.py - Tools for parsing CSV files
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

"""

import string, re, types
import csv

def ConvertDictionary( d , map = {}):
    """tries to convert values in a dictionary.

    if map contains 'default', a default conversion is enforced.
    For example, to force int for every column but column id,
    supply map = {'default' : "int", "id" : "str" }
    """

    rx_int = re.compile("^\s*[+-]*[0-9]+\s*$")
    rx_float = re.compile( "^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$" )

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

    for k,vv in d.items():

        if vv == None: continue
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
        except TypeError, msg:
            raise TypeError( "conversion in field: %s, %s" % (k,msg))

        try:
            if rx_int.match( v ):
                d[k] = int(v)
            elif rx_float.match( v ):
                d[k] = float(v)
        except TypeError, msg:
            raise TypeError( "expected string or buffer: offending value = '%s' " % str(v))
        except ValueError, msg:
            raise ValueError( "conversion error: %s, %s" % (msg, str(d)))
    return d

def GetMapColumn2Type( rows, ignore_empty = False, get_max_values = False ):
    """map fields to types based on rows.

    Preference is Int to Float to String.

    If get_max_values is set to true, the maximum values for integer
    columns are returned in a dictionary.
    """
    
    headers = rows[0].keys()
    map_column2type = {}
    is_full = {}

    max_values = {}
    
    for row in rows:
        for h in headers:

            if row[h] == "": continue
            
            is_full[h] = True

            if type(row[h]) == types.IntType:
                t = types.IntType
                if h not in max_values:
                    max_values[h] = int(row[h])
                else:
                    max_values[h] = max(h, int(row[h]))
                    
            elif type(row[h]) == types.FloatType:
                t = types.FloatType
            else:
                continue
            
            map_column2type[h] = t
            
    ignored = []
    for h in headers:
        if h not in map_column2type:
            if h in is_full or not ignore_empty:
                map_column2type[h] = types.StringType
            else:
                ignored.append( h )                

    if get_max_values:
        return map_column2type, ignored, max_values
    else:
        return map_column2type, ignored

##-----------------------------------------------------------
class CommentStripper:
    """iterator class for stripping comments from file.
    """
    
    def __init__(self, infile ):
        self.infile = infile

    def __iter__(self):
        return self
    
    def next(self):
        while 1:
            line = self.infile.next()
            if not line:
                raise StopIteration
            if line[0] != "#":
                return line

class DictReader( csv.DictReader ):
    """like csv.DictReader, but skip comments (lines starting with "#").
    """

    def __init__(self, infile, *args, **kwargs ):
        csv.DictReader.__init__( self,
                                 CommentStripper( infile ),
                                 *args, **kwargs )

#IMS - UnicodeDictReader to read Unicode encoded files
#Taken from http://stackoverflow.com/questions/1846135/python-csv-library-with-unicode-utf-8-support-that-just-works?rq=1
# WARNING: Does not strip # comment lines yet

class UnicodeCsvReader(object):
    def __init__(self, f, encoding="utf-8", **kwargs):
        self.csv_reader = csv.reader(f, **kwargs)
        self.encoding = encoding

    def __iter__(self):
        return self

    def next(self):
        # read and split the csv row into fields
        row = self.csv_reader.next() 
        # now decode
        return [unicode(cell, self.encoding) for cell in row]

    @property
    def line_num(self):
        return self.csv_reader.line_num

class UnicodeDictReader(csv.DictReader):
    def __init__(self, f, encoding="utf-8", fieldnames=None, **kwds):
        csv.DictReader.__init__(self, f, fieldnames=fieldnames, **kwds)
        self.reader = UnicodeCsvReader(f, encoding=encoding, **kwds)


class DictReaderLarge:
    """drop-in for csv.DictReader - handles very large fields
    
    Warning - minimal implementation - does not handle dialects
    """

    def __init__(self, infile, fieldnames, *args, **kwargs ):
        self.mFile = infile
        self.mFieldNames = fieldnames
        self.mNFields = len(fieldnames)

    def __iter__(self): 
        return self
    
    def next(self):
        
        line = self.mFile.next()
        if not line: raise StopIteration
        data = line[:-1].split("\t")
        assert len(data) == self.mNFields
        return dict( zip( self.mFieldNames, data ) )

######################################################################################
def ReadTable( lines,
               as_rows = True,
               with_header = True,
               ignore_incomplete = False,
               dialect = "excel-tab" ):
    """read a table from infile

    returns table as rows or as columns.
    If remove_incomplete, incomplete rows are simply ignored.
    """

    if type(lines) == types.FileType:
        lines = lines.readlines()

    lines = filter(lambda x: x[0] != "#", lines)

    if len(lines) == 0: return [], []

    if with_header:
        fields = lines[0][:-1].split("\t")
        del lines[0]
    else:
        fields = lines[0][:-1].split("\t")
        fields = map(str, range(len(fields)))

    nfields = len(fields)
    
    try:
        reader = csv.reader( lines.__iter__(), 
                             dialect = dialect )
    except TypeError:
        reader = csv.reader( lines.__iter__() )

    table = list(reader)
    
    if ignore_incomplete:
        table = [ x for x in table if len(x) == nfields ]
    else:
        for r, row in enumerate(table):
            if len(row) != nfields:
                if not ignore_incomplete:
                    raise ValueError( "missing elements in line %s, received=%s, expected=%s" % \
                                          (r, str(row),  str(fields)) )

                raise ValueError
        
    if not as_rows: table = zip( *table )
        
    return fields, table

######################################################################################
def ReadTables( infile, *args, **kwargs ):
    """read a set of csv tables. 

    Individual tables are separated by // on a single line.
    """
    
    lines = filter(lambda x: x[0] != "#", infile.readlines())
    chunks = filter( lambda x: lines[x][:2] == "//", range( len(lines )) )
    if not lines[-1].startswith("//"): chunks.append( len(lines) )

    class Result: pass
    result = []
    
    start = 0
    for end in chunks:
        
        fields, table = ReadTable( lines[start:end], *args, **kwargs )
        r = Result()
        r.mFields = fields
        r.mTable = table
        result.append(r)
        start = end + 1

    return result

######################################################################################
## group rows in table
def __DoGroup( rows, group_column, group_function, missing_value = "na"):

    values = []
    for x in range(len(rows[0])):
        if x == group_column:
            values.append( rows[0][x] )
        else:
            v = filter( lambda x: x != missing_value, map(lambda y: y[x], rows))
            if len(v) == 0:
                values.append( missing_value)
            else:
                values.append( group_function(map(lambda y: y[x], rows)))

    return values
    
def GroupTable( table,
                group_column = 0,
                group_function = min,
                missing_value = "na"):
    '''group table by *group_column*.

    The table need not be sorted.
    '''

    table.sort( lambda x,y: cmp(x[group_column], y[group_column]))

    rows = []
    last_value = None
    new_table = []
    
    for row in table:
        if row[group_column] != last_value:

            if last_value != None:
                new_table.append( __DoGroup( rows, group_column, group_function, missing_value ) )
                
            rows = []
            last_value = row[group_column]
            
        rows.append(row)

    if last_value != None:
        new_table.append( __DoGroup( rows, group_column, group_function, missing_value ) )

    return new_table

def getConvertedTable( table, columns, function = float,
                       skip_errors = False):
                   
    ## convert values to floats (except for group_column)
    ## Delete rows with unconvertable values
    new_table = []
    for row in table:
        skip = False
        for c in columns:
            try:
                row[c] = float(row[c])
            except ValueError, msg:
                if skip_errors:
                    skip = True
                    break
                else:
                    raise ValueError, msg
                
        if not skip: new_table.append(row)
        
    return new_table


