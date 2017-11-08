#!/usr/bin/env bash
#
# Script to find Python and R deps in this repository
#
# It uses snakefood to find the python dependencies
#

# exit when a command fails
#set -o errexit

# exit if any pipe commands fail
set -o pipefail

# exit when your script tries to use undeclared variables
#set -o nounset

# trace what gets executed
#set -o xtrace


# Global variables

SCRIPT_NAME="$0"
SCRIPT_OPTS="$@"
SCRIPT_FOLDER=$(dirname $0)
REPO_FOLDER=$(dirname ${SCRIPT_FOLDER})
TMP_SFOOD=$(mktemp)
TMP_GREP=$(mktemp)
TMP_MISC=$(mktemp)
TMP_DEPS=$(mktemp)

# dictionary to translate Python deps
declare -A PY_DEPS
PY_DEPS[Bio]="biopython"
PY_DEPS[MySQLdb]="mysqlclient"
PY_DEPS[SphinxReport]="ignore"
PY_DEPS[alignlib_lite]="alignlib-lite"
PY_DEPS[bs4]="beautifulsoup4"
PY_DEPS[bx]="bx-python"
PY_DEPS[configparser]="ignore"
PY_DEPS[drmaa]="python-drmaa"
PY_DEPS[future]="future"
PY_DEPS[ggplot]="ggplot"
PY_DEPS[jinja2]="jinja2"
PY_DEPS[lzo]="python-lzo"
PY_DEPS[matplotlib]="matplotlib"
PY_DEPS[metaphlan_utils]="ignore"
PY_DEPS[networkx]="networkx"
PY_DEPS[numpy]="numpy"
PY_DEPS[openpyxl]="openpyxl"
PY_DEPS[pandas]="pandas"
PY_DEPS[psycopg2]="psycopg2"
PY_DEPS[pyBigWig]="pybigwig"
PY_DEPS[pybedtools]="pybedtools"
PY_DEPS[pysam]="pysam"
PY_DEPS[pyximport]="ignore"
PY_DEPS[rdflib]="rdflib"
PY_DEPS[rpy2]="rpy2"
PY_DEPS[scipy]="scipy"
PY_DEPS[six]="six"
PY_DEPS[sklearn]="scikit-learn"
PY_DEPS[sqlalchemy]="sqlalchemy"
PY_DEPS[web]="web.py"
PY_DEPS[weblogolib]="python-weblogo"
PY_DEPS[yaml]="pyyaml"


# dictionary to translate R deps
declare -A R_DEPS
R_DEPS[BiSeq]="bioconductor-biseq"
R_DEPS[Biobase]="bioconductor-biobase"
R_DEPS[ChIPQC]="bioconductor-chipqc"
R_DEPS[DESeq2]="bioconductor-deseq2"
R_DEPS[DESeq]="bioconductor-deseq"
R_DEPS[DEXSeq]="bioconductor-dexseq"
R_DEPS[GMD]="r-gmd"
R_DEPS[HiddenMarkov]="r-hiddenmarkov"
R_DEPS[HilbertVis]="bioconductor-hilbertvis"
R_DEPS[Hmisc]="r-hmisc"
R_DEPS[IHW]="bioconductor-ihw"
R_DEPS[KEGGdb]="bioconductor-kegg.db"
R_DEPS[MASS]="r-mass"
R_DEPS[MEDIPS]="bioconductor-medips"
R_DEPS[RColorBrewer]="r-rcolorbrewer"
R_DEPS[RSQLite]="r-rsqlite"
R_DEPS[VennDiagram]="r-venndiagram"
R_DEPS[WGCNA]="r-wgcna"
R_DEPS[affy]="bioconductor-affy"
R_DEPS[amap]="r-amap"
R_DEPS[biomaRt]="bioconductor-biomart"
R_DEPS[coloc]="r-coloc"
R_DEPS[cummeRbund]="bioconductor-cummerbund"
R_DEPS[database]="ignore"
R_DEPS[dplyr]="r-dplyr"
R_DEPS[edgeR]="bioconductor-edger"
R_DEPS[flashClust]="r-flashclust"
R_DEPS[gcrma]="bioconductor-gcrma"
R_DEPS[genome_file]="ignore"
R_DEPS[ggplot2]="r-ggplot2"
R_DEPS[gplots]="r-gplots"
R_DEPS[gridExtra]="r-gridextra"
R_DEPS[grid]="r-gridbase"
R_DEPS[gtools]="r-gtools"
R_DEPS[hpar]="bioconductor-hpar"
R_DEPS[limma]="bioconductor-limma"
R_DEPS[maSigPro]="bioconductor-masigpro"
R_DEPS[mapdata]="r-mapdata"
R_DEPS[maps]="r-maps"
R_DEPS[metagenomeSeq]="bioconductor-metagenomeseq"
R_DEPS[optparse]="r-optparse"
R_DEPS[plotrix]="r-plotrix"
R_DEPS[plyr]="r-plyr"
R_DEPS[pvclust]="r-pvclust"
R_DEPS[qqman]="r-qqman"
R_DEPS[reshape2]="r-reshape2"
R_DEPS[reshape]="r-reshape"
R_DEPS[rtracklayer]="bioconductor-rtracklayer"
R_DEPS[samr]="r-samr"
R_DEPS[scales]="r-scales"
R_DEPS[sciplot]="r-sciplot"
R_DEPS[siggenes]="bioconductor-siggenes"
R_DEPS[simpleaffy]="bioconductor-simpleaffy"
R_DEPS[sleuth]="r-sleuth"
R_DEPS[snow]="r-snow"
R_DEPS[spp]="ignore"
R_DEPS[vegan]="r-vegan"
R_DEPS[wasabi]="r-wasabi"
R_DEPS[zinba]="ignore"


# dictionary to translate Misc deps
declare -A MISC_DEPS
MISC_DEPS[ascp]="ignore"
MISC_DEPS[bamToBed]="bedtools"
MISC_DEPS[bedGraphToBigWig]="ucsc-bedgraphtobigwig"
MISC_DEPS[bedToBigBed]="ucsc-bedtobigbed"
MISC_DEPS[bedtools]="bedtools"
MISC_DEPS[bits_test]="ignore"
MISC_DEPS[cat]="coreutils"
MISC_DEPS[gat-run.py]="gat"
MISC_DEPS[gdc-client]="gdc-client"
MISC_DEPS[grep]="grep"
MISC_DEPS[gzip]="ignore" # gzip not available in conda
MISC_DEPS[intersectBed]="bedtools"
MISC_DEPS[java]="ignore"
MISC_DEPS[paste]="coreutils"
MISC_DEPS[prefetch]="ignore"
MISC_DEPS[python]="ignore"
MISC_DEPS[rename]="util-linux"
MISC_DEPS[rm]="coreutils"
MISC_DEPS[samtools]="samtools"
MISC_DEPS[sort]="coreutils"
MISC_DEPS[srapath]="sra-tools"
MISC_DEPS[tar]="tar"
MISC_DEPS[tr]="coreutils"
MISC_DEPS[wget]="wget"
MISC_DEPS[wigToBigWig]="ucsc-wigtobigwig"
MISC_DEPS[zcat]="ignore" # gzip not available in conda


# function to report issues and exit
report_problem() {
   echo
   echo $1
   echo
   echo " Aborting. "
   exit 1
}


# function to find python imports
# output will go to TMP_SFOOD
find_python_imports() {

sfood $1 2>&1 \
 | grep 'WARNING     :   ' \
 | grep -v Line \
 | awk '{print $3;}' \
 | grep -v '^_.*' \
 | sed 's/\..*//g' \
 | egrep -v 'CGAT|XGram|builtins|corebio|pylab|xml' \
 | sort -u \
 >> ${TMP_SFOOD}

}


# function to find R imports
# output will go to TMP_GREP
find_r_imports() {

grep -i 'library(' -r $1 \
 | grep -v Binary \
 | sed -e 's/\(.*\)library\(.*\)/\2/' \
 | sed 's/[()"&,.%'\'']//g' \
 | sed 's/\\n$//g' \
 | egrep '^[a-zA-Z]{2,}' \
 | sort -u \
 >> ${TMP_GREP}

}


# function to find misc programs
# output will go to TMP_MISC
find_misc_programs() {

   # Not sure why, but this is required to run properly
   FIND_DIR=$1

   # will use specific py3 env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate py3-basic

   TMP_EXT=$(mktemp)
   find ${FIND_DIR} -iname "*.py" > ${TMP_EXT}

   for code in `cat ${TMP_EXT}` ;
   do

      python ${REPO_FOLDER}/scripts/cgat_check_deps.py -s ${code} >> ${TMP_MISC}

   done

   # return unique names
   awk '/^Program/ {print $2}' ${TMP_MISC} | sort -u > ${TMP_EXT}
   cp ${TMP_EXT} ${TMP_MISC}

   # revert to original env
   source deactivate
   source /ifs/apps/conda-envs/bin/activate snakefood

   # clean up tmp file
   rm ${TMP_EXT}

}


# function to display help message
help_message() {
   echo
   echo " Scans this repository to look for conda dependencies."
   echo
   echo " To get the dependencies for all scripts, run:"
   echo " ./cgat_conda_deps.sh --all"
   echo
   echo " To get the dependencies for production scripts only, run:"
   echo " ./cgat_conda_deps.sh --production"
   echo
   exit 1
}


# the script starts here

if [[ $# -eq 0 ]] ; then

   help_message

fi

# variable to choose scope
ALL=1

# parse input parameters

while [[ $# -gt 0 ]]
do
key="$1"

case $key in

    --help)
    help_message
    ;;

    --all)
    ALL=1
    shift
    ;;

    --production)
    ALL=0
    shift
    ;;

    *)
    help_message
    ;;

esac
done

# requirement: snakefood
source /ifs/apps/conda-envs/bin/activate snakefood

# initialize temp files
cat /dev/null > ${TMP_SFOOD}
cat /dev/null > ${TMP_GREP}
cat /dev/null > ${TMP_MISC}
cat /dev/null > ${TMP_DEPS}

# find deps depending on given input
if [[ ${ALL} -eq 1 ]] ; then

   # Python
   find_python_imports "${REPO_FOLDER}/CGAT"

   # R
   find_r_imports "${REPO_FOLDER}/CGAT"
   find_r_imports "${REPO_FOLDER}/R"

   # Misc
   find_misc_programs "${REPO_FOLDER}/CGAT"

else
   # Use tmp folder
   TMP_D=$(mktemp -d)/cgat
   mkdir -p ${TMP_D}

   # Bring production scripts to tmp folder
   grep CGAT ${REPO_FOLDER}/MANIFEST.in \
    | egrep -v '#|install|exclude|h$|c$|cpp$|pxd$|pyx' \
    | awk '{print $2;}' \
    > ${TMP_D}/grep-patterns

   for f in `find ${REPO_FOLDER}/CGAT | grep -f ${TMP_D}/grep-patterns` ; do cp $f ${TMP_D}/ ; done
   find ${REPO_FOLDER}/CGAT -regex '.*pyx.*' -exec cp {} ${TMP_D}/ \;

   # Also bring associated module files
   for f in `find ${TMP_D}` ; do awk '/^import CGAT/ {print $2}' $f 2>/dev/null ; done \
    | sort -u \
    | egrep -v 'scripts|NCL|GeneModelAnalysis' \
    | sed 's/\./\//g' \
    > ${TMP_D}/imported-modules

   for m in `cat ${TMP_D}/imported-modules` ; do cp ${REPO_FOLDER}/$m.py ${TMP_D}/ ; done

   # Python
   find_python_imports "${TMP_D}"

   # R
   find_r_imports "${TMP_D}"

   # Misc
   find_misc_programs "${TMP_D}"

   # clean up
   [[ -n "${TMP_D}" ]] && rm -rf "${TMP_D}"

fi


### create header of env file ###

echo
echo "# output generated by ${SCRIPT_NAME} ${SCRIPT_OPTS}"
echo "# on `date`"
echo
echo "name: cgat-s"
echo
echo "channels:"
echo "- bioconda"
echo "- conda-forge"
echo "- defaults"
echo
echo "dependencies:"


### process python deps ###

for pkg in `cat ${TMP_SFOOD}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${PY_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${PY_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${PY_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "? "$pkg >> ${TMP_DEPS}
   fi

done

# print Python section
echo "# python dependencies"
echo "- python"

# Add others manually:
echo "- cython" >> ${TMP_DEPS}
[[ ${ALL} -eq 1 ]] && echo "- nose" >> ${TMP_DEPS}
[[ ${ALL} -eq 1 ]] && echo "- pep8" >> ${TMP_DEPS}
echo "- setuptools" >> ${TMP_DEPS}
echo "- pyyaml" >> ${TMP_DEPS}

# Print them all sorted
sed 's/^- bx-python/# WARNING: bx-python is Py2 only but "pip install bx-python" works with Py3/g' ${TMP_DEPS} \
 | sort -u


### process R deps ###

cat /dev/null > ${TMP_DEPS}

for pkg in `cat ${TMP_GREP}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${R_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${R_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${R_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "- "$pkg >> ${TMP_DEPS}
   fi

done

# print R section
R_DEPS_SIZE=$(stat --printf="%s" ${TMP_DEPS})

if [[ ${R_DEPS_SIZE} -gt 0 ]] ; then

   echo "# R dependencies"
   echo "- r-base"
   sort -u ${TMP_DEPS} | grep '^\- r'

   # show bioconductor deps
   sort -u ${TMP_DEPS} | grep '^\- bioconductor'

   # show missing deps
   sort -u ${TMP_DEPS} | grep '^? '

fi


### process misc programs ###

cat /dev/null > ${TMP_DEPS}

for pkg in `cat ${TMP_MISC}` ;
do

   # Reference:
   # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
   if [[ ${MISC_DEPS[${pkg}]+_} ]] ; then
      # found
      [[ "${MISC_DEPS[${pkg}]}" != "ignore" ]] && echo "- "${MISC_DEPS[${pkg}]} >> ${TMP_DEPS}
   else
      # not found
      echo "? "$pkg >> ${TMP_DEPS}
   fi

done

# print misc section

# Add these manually, as they are required but don't use the 'statement' variable to run them
echo "- nomkl" >> ${TMP_DEPS}
echo "- gcc" >> ${TMP_DEPS}
echo "- zlib" >> ${TMP_DEPS}
echo "- libpng" >> ${TMP_DEPS}

echo "# Misc dependencies"
egrep -v 'CGATparameter|perm|--' ${TMP_DEPS} \
 | sort -u


# Remove temp files
rm ${TMP_SFOOD}
rm ${TMP_GREP}
rm ${TMP_MISC}
rm ${TMP_DEPS}

