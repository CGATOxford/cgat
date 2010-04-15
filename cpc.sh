#!/bin/bash

######################################################################
# 
# Arguments
# 
arg_input_seq=$1

arg_output_result=$2

arg_working_dir=$3

arg_output_evd_plot_feat_base=$4

arg_blast_db=$5

######################################################################
# 
# Constants: for the directory settings
# 

MYOWN_LOC=/net/cpp-group/tools/cpc-0.9/bin;		# XXX: some hacking:)
DATA_DIR="$MYOWN_LOC/../data"
LIB_DIR="$MYOWN_LOC/../libs"

## location of the farm script to send jobs to the cluster
APP_FARM="/net/cpp-group/scripts/farm.py --split-at-regex=^>(\S+) --chunksize=50 --log=${arg_working_dir}/farm.log --cluster-queue=medium_jobs.q --cluster-priority=-10 --tmpdir=${arg_working_dir}"

if [ -n "$arg_blast_db" ]; then
    m_blast_db=$arg_blast_db;
else 
    m_blast_db=$DATA_DIR/prot_db;
fi 

echo "using blast database ${m_blast_db}"

m_framefinder_model=$DATA_DIR/framefinder.model

m_libsvm_model0=$DATA_DIR/libsvm.model0 # standard
m_libsvm_model=$DATA_DIR/libsvm.model # Prob
m_libsvm_model2=$DATA_DIR/libsvm.model2	# Prob + weighted version
m_libsvm_range=$DATA_DIR/libsvm.range

c_extract_blast_feat="$MYOWN_LOC/extract_blastx_features.pl"
c_extract_ff_feat="$MYOWN_LOC/extract_framefinder_feats.pl"
c_add_missing_entries="$MYOWN_LOC/add_missing_entries.pl"
c_feat2libsvm="$MYOWN_LOC/feat2libsvm.pl"
c_lsv_cbind="$MYOWN_LOC/lsv_cbind.pl"
c_join_columns="$MYOWN_LOC/join_columns.pl"
c_predict="$MYOWN_LOC/predict.pl"
c_generate_plot_feats="$MYOWN_LOC/generate_plot_features.pl"
c_split_plot_feats="$MYOWN_LOC/split_plot_features_by_type.pl"
c_index_blast_report="$MYOWN_LOC/make_blast_report_index.pl"

# 
# Constants: for the remote blast
# 
REMOTE_BLAST_HOST="162.105.250.200" # New LangChao
REMOTE_BLAST_MIN_SIZE=4000	# 4k
c_blast_smp_client="$MYOWN_LOC/server/client.pl"

if test ! -d "${arg_working_dir}"; then
    mkdir -p ${arg_working_dir}; \
    echo "creating working dir ${arg_working_dir}"; \
fi


############################################################
# 
# Step 0: detect necessary applications
# 
APP_BLAST=`which blastall 2> /dev/null`
test -x "$APP_BLAST" || (echo "Can't find blastall on your path, check it!" > /dev/stderr && exit 1)
## reset to blastall found in path, this will map it for the nodes
APP_BLAST="blastall"

APP_FF=`which framefinder 2> /dev/null`	# FF == FrameFinder
if test ! -x "$APP_FF"; then
    APP_FF=$LIB_DIR/framefinder;
    if test ! -x "$APP_FF"; then
	APP_FF=$LIB_DIR/estate/bin/framefinder;
	if test ! -x "$APP_FF"; then
	    echo "Can't find framefinder on your path or my own directory, quitting..." > /dev/stderr
	    exit 1
	fi;
    fi
fi

APP_SVM_SCALE="$LIB_DIR/libsvm/libsvm-2.81/svm-scale"
test -x "$APP_SVM_SCALE" || (echo "Can't find svm-scale on your path, eheck it!" > /dev/stderr && exit 1)

APP_SVM_PREDICT="$LIB_DIR/libsvm/libsvm-2.81/svm-predict"
APP_SVM_PREDICT2="$LIB_DIR/libsvm/libsvm-2.81/svm-predict2"
test -x "$APP_SVM_PREDICT" || (echo "Can't find svm-predict on your path, eheck it!" > /dev/stderr && exit 1)
test -x "$APP_SVM_PREDICT2" || (echo "Can't find svm-predict2 on your path, eheck it!" > /dev/stderr && exit 1)

APP_BLAST2TAB=$LIB_DIR/blast2table.pl
if test ! -x "$APP_FF"; then
    chmod +x $APP_BLAST2TAB 2> /dev/null;
    if test ! -x "$APP_FF"; then
	echo "Can't find blast2table.pl on my own package, something goes wrong..." > /dev/stderr
	exit 1
    fi
fi

# detect the BLAST database, failsafe
if test ! -e "${m_blast_db}.phr" && test ! -e "${m_blast_db}.pal"; then
    echo "Can't find protein db ${m_blast_db}.{phr,pal} under $DATA_DIR, pls check..." > /dev/stderr
    exit 1
fi


# Step 1: run blastx & framefinder

# BLASTX settings: Combining the BLAST and Frith2006(PLoS & RNA) protocols
# XXX: the remote server will NOT use their own settings...
blast_opts="-S 1";              # only the same strand
blast_opts="$blast_opts -e 1e-10"; # as a quick setting (BLAST 9.3.2)
blast_opts="$blast_opts -g F";  # un-gapped blast (Frith2006, PLoS)
blast_opts="$blast_opts -f 14"; # Neighborhood word threshold score, default=12 (BLAST 9.3.2)
blast_opts="$blast_opts -a 2";  # 2 CPUs, boost the performance

blast_opts="$blast_opts -d $m_blast_db"	# database settings

# Framefinder settings
ff_opts="-r False -w $m_framefinder_model /dev/stdin"

# Entry the working space...
old_pwd=`pwd`

# cd $arg_working_dir || (echo "Can't enter the working space ($arg_working_dir), quitting...." > /dev/stderr && exit 2)

# Determine the right mode (local or remote) for running BLAST
input_seq_size=`stat -Lc "%s" $arg_input_seq`;

# local version
(cat $arg_input_seq | $APP_FARM $APP_BLAST -p blastx $blast_opts | tee $arg_working_dir/blastx.bls | perl $APP_BLAST2TAB | tee $arg_working_dir/blastx.table | perl $c_extract_blast_feat ) > $arg_working_dir/blastx.feat1 &

(cat $arg_input_seq | $APP_FF $ff_opts | tee $arg_working_dir/ff.fa1 | perl $c_extract_ff_feat ) > $arg_working_dir/ff.feat &

wait;

# a quick fix: adding possible missing entries due to blastx
cat $arg_input_seq | perl $c_add_missing_entries $arg_working_dir/blastx.feat1 > $arg_working_dir/blastx.feat
# Quick fix: remove redunancy \r in the ff.fa
cat $arg_working_dir/ff.fa1 | tr -d '\r' > $arg_working_dir/ff.fa

############################################################
# 
# Step 2: prepare data for libsvm
# 

# 1       2             3              4         5           6
# QueryID hit_seq_count hit_HSP_count  hit_score frame_score frame_score2

perl $c_feat2libsvm -c 2,4,6 NA NA $arg_working_dir/blastx.feat > $arg_working_dir/blastx.lsv &

# 1       2             3              4         5
# QueryID CDSLength     Score          Used     Strict
perl $c_feat2libsvm -c 2,3,4,5 NA NA $arg_working_dir/ff.feat > $arg_working_dir/ff.lsv &
wait;

perl -w $c_lsv_cbind $arg_working_dir/blastx.lsv $arg_working_dir/ff.lsv > $arg_working_dir/test.lsv
$APP_SVM_SCALE -r $m_libsvm_range $arg_working_dir/test.lsv > $arg_working_dir/test.lsv.scaled

############################################################
# 
# Step 3: do prediction
# 

$APP_SVM_PREDICT2 $arg_working_dir/test.lsv.scaled $m_libsvm_model0 $arg_working_dir/test.svm0.predict > $arg_working_dir/test.svm0.stdout 2> $arg_working_dir/test.svm0.stderr

cat $arg_working_dir/test.svm0.predict  | perl -w $c_predict $arg_input_seq > $arg_output_result

############################################################
# 
# Step 4: generate the output features for web-visualization
# 

output_plot_feat_homo=${arg_output_evd_plot_feat_base}.homo
output_plot_feat_orf=${arg_output_evd_plot_feat_base}.orf

cat $arg_working_dir/blastx.feat | perl -w $c_generate_plot_feats $arg_working_dir/blastx.table $arg_working_dir/ff.fa | perl -w $c_split_plot_feats $output_plot_feat_homo $output_plot_feat_orf &

perl -w $c_index_blast_report $arg_working_dir/blastx.bls > $arg_working_dir/blastx.index &

wait;

############################################################
# 
# Step 5: make some clean-up...
# 

rm -rf $arg_working_dir/blastx.feat1
rm -rf $arg_working_dir/ff.fa1

# leaving out the working space...
cd $old_pwd
