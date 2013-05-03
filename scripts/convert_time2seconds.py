"""convert output of the time command into seconds:

input: for example something like

1h33m47.32s

"""

import os, sys, string, re

if __name__ == "__main__":

    for line in sys.stdin:

        if line[0] == "#": continue

        data = string.split(line[:-1], "\t")

        new_fields = []

        for field in data:

            x = re.match("^(\d+)h(\d+)m([\d.]+)s", field)
            if not x:
                x = re.match("^(\d+)m([\d.]+)s", field)
                if x:
                    minutes, seconds = x.groups()
                else:
                    new_fields.append(field)
                    continue
                hours = "0"
            else:
                hours, minutes, seconds = x.groups()


            new_fields.append( "%.2f" % (string.atof(hours) * 3600.0 + string.atof(minutes) * 60.0 + string.atof(seconds)))

        print string.join(new_fields, "\t")


    



    

