import sys, re, string, os, getopt, time

USAGE = """python %s < stdin > stdout

substitute tokens

OPTIONS:
-c, --create            create substitution list
-a, --apply             apply substitution list
-i, --invert            pairs of substitution patterns is to,from
-r, --regex-token       regular expression for tokens (has to create one pair of brackets)
-p, --pattern-sub       pattern for substitution
-m, --multiple          do multiple substitutions per row
-o, --columns-token     substitute tokens in columns
-e, --echo              echo susbstituted column
-f, --filter            remove lines not matching
-y, --reverse-filter    remove lines matching (filters for errors)
-n, --inplace           do inplace subsitutions of all files on command line
-b, --backup            keep backup (with ending .bak)
-s, --select-rows       regular expression for rows to use.
-x, --extended          replace not just with second column in map, but all columns.
-k, --keep              keep column that is substituted
--keep-header           do not apply transformation to header
""" % sys.argv[0]

param_long_options = ["verbose=",
                      "create=", "regex-token=", "pattern-sub=",
                      "apply=", "invert", "multiple", "columns-token=", "echo",
                      "filter", "inplace", "backup", "select-rows=", "extended", "keep",
                      "reverse-filter", "keep-header" ]
param_short_options = "c:r:p:a:imo:fnbs:v:x"

param_create = None
param_regex_token = None
param_pattern_sub = "%s"
param_apply = None
param_invert = None
param_multiple = None
param_columns_token = None

param_filter = False
param_reverse_filter = False

param_inplace = False
param_backup = False

param_loglevel = 1

param_regex_rows=None

##replacement options
param_echo = False
param_extended = False
param_keep = False

param_keep_header = False

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)
                                      

    except getopt.error, msg:
        print USAGE, msg
        sys.exit(1)

    for o,a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("-c", "--create"):
            param_create = a
        elif o in ("-r", "--regex-token"):
            param_regex_token = re.compile( a )
        elif o in ("-p", "--pattern-sub"):
            param_pattern_sub = a
        elif o in ("-a", "--apply"):
            param_apply = a
        elif o in ("-x", "--extended"):
            param_extended = True
        elif o in ("-i", "--invert"):
            param_invert = 1
        elif o in ("-m", "--multiple"):
            param_multiple = 1
        elif o in ("-o", "--columns-token"):
            if a != "all":
                param_columns_token = map(lambda x: int(x) - 1,string.split(a, ","))
            else:
                param_columns_token = a
        elif o in ("-e", "--echo"):
            param_echo = True
        elif o in ("-k", "--keep"):
            param_keep = True
        elif o in ("-f", "--filter"):
            param_filter = True
        elif o in ("-y", "--reverse-filter"):
            param_reverse_filter = True
        elif o in ("-n", "--inplace"):
            param_inplace = True
        elif o in ("-b", "--backup"):
            param_backup = True
        elif o in ("-s", "--select-rows"):
            param_regex_rows = re.compile(a)
        elif o in ("--keep-header"):
            param_keep_header = True

    if param_inplace and len(args) == 0:
        raise "please supply file name(s) for inplace substitution."

    file_id = 0

    keys = {}

    if param_apply:
        infile = open(param_apply, "r")
        for line in infile:
            if line[0] == "#": continue

            d = line[:-1].split("\t")
            try:
                a, b = d[:2]
            except ValueError:
                print "# invalid map skipped in line: %s" % line
                continue
            
            if param_invert:            
                a, b = b, a
                if param_extended:
                    b = "\t".join(d[0] + d[2:])
            else:
                if param_extended:
                    b = "\t".join(d[1:])
                    
            if not keys.has_key( a ):
                keys[a] = []

            if param_keep:
                b = a + "\t" + b

            keys[a].append( b )
            
    files = args
    
    if not param_inplace and len(args) == 0:
        files = ["-"]
        
    for file in files:

        close_infile = False
        close_outfile =False
        if file == "-":
            infile = sys.stdin
            outfile = sys.stdout
        else:
            if param_inplace:
                os.rename( file, file + ".bak" )
                infile = open( file + ".bak", "r")
                outfile = open( file, "w")
                close_infile = True
                close_outfile = True
            else:
                infile = open( file, "r")
                outfile = sys.stdout
                close_infile = True

        first = True
        
        for line in infile:
            if line[0] == "#":
                outfile.write(line)
                continue

            if first:
                first = False
                if param_keep_header:
                    outfile.write(line)
                    continue

            if param_regex_rows:
                if param_regex_rows.search( line ):
                    outfile.write(line)
                    continue
                
            new_lines = []
            if param_regex_token:
                r = param_regex_token.search( line[:-1] )
                while r:
                    key = r.group(1)
                    if key not in keys:
                        if param_create:
                            keys[key] = [param_pattern_sub % str(len(keys))]
                        else:
                            new_lines.append(line[:-1])
                            break

                    for k in keys[key]:
                        new_lines.append( line[:r.start(1)] + k + line[r.end(1):-1] )

                    if param_multiple:
                        r = param_regex_token.search( line[r.end(1):-1] )
                    else:
                        break
                else:
                    if not param_filter:
                        new_lines.append(line[:-1])

            elif param_columns_token:
                data = line[:-1].split("\t")
                if param_columns_token == "all":
                    columns = range( len(data))
                else:
                    columns = param_columns_token
                keep = not param_reverse_filter
                first_multiple = True
                for c in columns:
                    k = data[c]
                    if k in keys:
                        if len(keys[k]) > 1:
                            if not first_multiple: 
                                raise "warning: could not substitute multiple keys for %s in multiple columns in line: %s" % (k, line)
                            first_multiple = False
                        for v in keys[k]:
                            if param_echo: data.append(data[c])
                            ## multiple substitutions: write data now
                            data[c] = v
                            if keep:
                                new_lines.append( string.join( data, "\t" ) )
                        keep = False
                    else:
                        if param_create:
                            keys[k] = [param_pattern_sub % str(len(keys))]
                            data[c] = keys[k][0]
                        elif param_filter:
                            keep = False
                        elif param_reverse_filter:
                            keep = True
                if keep:
                    new_lines.append(string.join(data, "\t"))

            elif param_apply:
                for key in keys:
                    for k in keys[key]:
                        line = line.replace( key, k )
                new_lines.append( line[:-1] )

            if new_lines:
                outfile.write(string.join(new_lines, "\n") + "\n")

        if param_create:
            create_file = open(param_create, "w")
            for key in keys:
                for k in keys[key]:
                    create_file.write("%s\t%s\n" % (key, str(k)))
            create_file.close()

        if close_outfile: outfile.close()
        if close_infile: infile.close()
    
        if param_inplace and not param_backup:
            os.remove( file + ".bak" )
