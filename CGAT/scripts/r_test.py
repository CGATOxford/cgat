'''
r_test.py - t-test and wilcoxon tests
=====================================

:Tags: Python

Purpose
-------

Interface to R to run tests on a vector of values.

Usage
-----

Example::

   python r_test.py --help

Type::

   python r_test.py --help

for command line help.

Command line options
--------------------

'''
import sys
from rpy2.robjects import r as R
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: r_test.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      help="method to use [t-test=t-test,wilcox=wilcox]",
                      choices=("t-test", "wilcox"))
    parser.add_option("-1", "--infile", dest="filename_input", type="string",
                      help="input filename with vector of values.")
    parser.add_option("-2", "--infile2", dest="filename_input2", type="string",
                      help="input filename with vector of values.")
    parser.add_option("--header-names", dest="header", type="string",
                      help="""header of value column [default=%default].""")

    parser.set_defaults(
        method="t-test",
        filename_input=None,
        header="value",
    )

    (options, args) = E.Start(parser,
                              add_pipe_options=True)

    if options.filename_input:
        infile = IOTools.openFile(options.filename_input, "r")
    else:
        infile = sys.stdin

    values, errors = IOTools.ReadList(infile,
                                      map_function=float)
    if options.filename_input:
        infile.close()

    if errors:
        E.warn("errors in input: %s" % ";".join(map(str, errors)))

    kwargs = {}
    xargs = []
    for arg in args:
        if "=" in arg:
            key, value = arg.split("=")
            kwargs[key] = value
        else:
            xargs.append(arg)

    if options.filename_input2:
        infile = IOTools.openFile(options.filename_input2, "r")
        values2, errors2 = IOTools.ReadList(infile,
                                            map_function=float)
        infile.close()
    else:
        values2 = None

    stat = Stats.Summary(values)

    power, diff_at_power95 = None, None
    if options.method == "t-test":
        if values2:
            result = R.t_test(values, values2, *xargs, **kwargs)
        else:
            result = R.t_test(values, *xargs, **kwargs)
            # compute power of test
            power = R.power_t_test(n=len(values),
                                   delta=abs(stat["mean"]),
                                   sd=stat["stddev"],
                                   sig_level=0.05)['power']
            diff_at_power95 = R.power_t_test(n=len(values),
                                             power=0.95,
                                             sd=stat["stddev"],
                                             sig_level=0.05)['delta']

    if options.method == "wilcox":
        result = R.wilcox_test(values, *xargs, **kwargs)

    options.stdout.write("%s\t%s\n" % ("key", options.header))

    for key, value in sorted(result.items()):
        if key == "data.name":
            continue
        if key == "p.value":
            options.stdout.write("%s\t%5.2e\n" % (str(key), value))
        else:
            options.stdout.write("%s\t%s\n" % (str(key), str(value)))

    for key, value in list(stat.items()):
        options.stdout.write("%s\t%s\n" % (str(key), str(value)))

    if power:
        options.stdout.write("1-power\t%5.2e\n" % (1.0 - power))
        options.stdout.write("diff_at_power95\t%f\n" % diff_at_power95)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
