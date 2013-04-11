# CGAT, Steve, March 2013 with credit to Tim for implementation suggestions.
#
# This script allows a specificied function from a specified python module 
# to be executed on the cluster with specified parameters.
#
# A statement is specified in the normal way i.e.:
#
# statement = ''' python %(scriptsdir)s/run_function.py 
#                        -p infile,outfile,additional_param1
#                        -m modulefile 
#                        -f function '''
#
# P.run()
#
# If the module is in your $PYTHONPATH you can just name it 
# directly. i.e "Pipeline" would specifiy Pipeline.py
#

from optparse import OptionParser
import sys, os

def main(argv = None):
    
    # Parse the options
    parser = OptionParser()   
    parser.add_option("-p", "--params", dest="params", type="string",
                      help="comma separated list of addtional parameter strings")
    parser.add_option("-m", "--module", dest="module", type="string",
                      help="the full path to the module file", default = None)
    parser.add_option("-f", "--function", dest="function", type="string",
                      help="the module function", default = None)
    (options, args) = parser.parse_args()

    # Check a module and function have been specified
    if not options.module or not options.function:
        raise ValueError("Both a function and Module must be specified")

    # If a full path was given, add this path to the system path
    location = os.path.dirname(options.module)
    if location != "":
        sys.path.append(location)
   
    # Establish the module name, accomodating cases where the
    # .py extension has been included in the module name
    module_name = os.path.basename(options.module)
    if module_name.endswith(".py"):
        module_base_name = module_name[:-3]
    else:
        module_base_name = module_name

    # Import the specified module and map the specified fuction
    module = __import__(module_base_name)
    function = getattr(module, options.function)

    # Parse the parameters into an array
    params = [param.strip() for param in options.params.split(",")]
              
    # Make the function call
    function(params)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
