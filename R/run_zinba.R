
# Parse all Rscript arguments
parse_command_args <- function()
{
		l <- list()
	 	a <- commandArgs(trailingOnly=T)
	 	for (arg in a)
	 	{
			bits <- unlist(strsplit(arg,"="))
			name=bits[1]
			value=bits[2]
			l[[name]] <- value
	 	}
	  return(l)
}

cargs <- parse_command_args()

# Check all required arguments are present
if(!(sum(names(cargs)=="outname") & sum(names(cargs)=="fdr") & sum(names(cargs)=="chip") & sum(names(cargs)=="input")       ))
{
	print("chip, input, outname and fdr are required")
	q()
} else {
  outname = as.character(cargs["outname"])
  chipfile = as.character(cargs["chip"])
  inputfile = as.character(cargs["input"])
  fdr = as.numeric(cargs["fdr"])
  print(paste("outname:",outname,"chip:",chipfile,"input",inputfile,"fdr",fdr))
}

#load the library
library('zinba')

generateAlignability(
   mapdir=,      #mappability directory from unpacked mappability files
   outdir=,      #directory for processed files, used later in analysis
   athresh=,     #number of hits per read allowed during mapping process
   extension=,   #average fragment library length
   twoBitFile=,  #path to downloaded genome build file in .2bit format
)

basealigncount(
   inputfile=,  #mapped sample reads
   outputfile=, # output path
   extension=,  #average fragment library length 
   filetype=,   #either "bed", "bowtie", or "tagAlign"
   twoBitFile=, #path to downloaded genome build file in .2bit format
)

zinba(
  refinepeaks=,   #refine peaks?  1 for yes, 0 for no
  seq=,           #path to mapped experimental reads
  input=,         #path to mapped input reads if available (default is "none")
  filetype=,      #either 'bed', 'bowtie', or 'tagAlign'
  threshold=,     #FDR threshold, default is 0.05 
  align=,         #path to alignability directory
  numProc=,       #number of CPUs to use, must be less than max available (default 1)
  twoBit=,        #path to genome build in .2bit format
  outfile=,       #prefix for outputted files
  extension=,     #average fragment library length (size selected)
  basecountfile=, #path to basecount file if refinepeaks is 1
  broad=,         #broad setting, TRUE or FALSE (default)
  printFullOut=,  #print original data with enrichment estimates, 1  for yes (more space required), 0 for no (default)
  interaction=,   #whether or not to considering interaction during model selection, TRUE (default) or FALSE
  mode=,          #either "peaks" for peak calling (default) or "CNV" for calling likely amplified CNV regions for reads in "seq" (input reads are best)
  FDR=            #either TRUE (default) or FALSE.  If false, then uses posterior probability to threshold peaks using 1-threshold 
)


