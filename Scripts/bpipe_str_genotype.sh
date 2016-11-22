#!/bin/bash
#
TOOLS_DIR='/path/to/tools(bpipe)'
bpipe_ver='0.9.8.7'
BPIPE_EXE="${TOOLS_DIR}/bpipe-${bpipe_ver}/bin/bpipe"

STRSEQ_DIR='/path/to/strseq'
SCRIPT_DIR="${STRSEQ_DIR}/Scripts"
PIPES_DIR="${STRSEQ_DIR}/Scripts/bpipes"

# Function to print help
print_usage()
{
	   echo "Usage: $(basename "$0") -p {OSSeq_pool} [-r {rd_len} -k {flank_size} -v {snv_filter} -f {file set} ]";
	   echo "Where -p pool file must be specified (eg OS0037)";
	   echo "      -r read length, default = 150";
	   echo "      -k flank size, default = 15";
	   echo "      -v R2 snv filter [all|alt], default = alt";     
	   echo "      -f override file(s) to be run in pipeline, default = '*.ptag.st.bam'";
	   echo "      -d debug (do not delete intermediate files)";
	   return
}

# set default parameters
fset=""; rdlen=""; flank=""; snvs="";
default_fset='*.ptag.st.bam'
default_rdlen=150
default_flank=15
default_snvs='alt'
clean='Y'
re_int='^[0-9]+$'

# Parse command line options
OPTIND=1
while getopts "p:r:k:v:f:dh" OPT
do
  case "$OPT" in
    p) pool="$OPTARG";;
    r) rdlen="$OPTARG";;
    k) flank="$OPTARG";;
    v) snvs="$OPTARG";;
    f) fset="$OPTARG";;
    d) clean="N";;
    h) print_usage; exit 1;;
   \?) print_usage; exit 1;;
    :) echo "Option -$OPTARG requires an argument."; print_usage; exit 1;;
  esac
done

if [ "$pool" == "" ]; then
  echo "OSSeq pool not specified"
  print_usage; 
  exit 2
fi

if [ "$rdlen" == "" ]; then
  rdlen=$default_rdlen
elif ! [[ $rdlen =~ $re_int ]]; then
  echo "Read length provided ($rdlen) is not an integer"; exit 1
fi

if [ "$flank" == "" ]; then
  flank=$default_flank
elif ! [[ $flank =~ $re_int ]]; then
  echo "Flank size provided ($flank) is not an integer"; exit 1
fi

if [ "$snvs" == "" ]; then
  snvs=$default_snvs
fi

if [ "$fset" == "" ]; then
  fset=$default_fset
fi

export FLANK_SIZE=$flank  

# Turn on extended globbing (allowing regex in file matching patterns)
shopt -s extglob

# Define stages to execute
ext=${fset: -3}
if [ "$ext" != 'bam' ]; then
  echo "Input file set must have .bam extension"
  print_usage;
  exit 2
fi

echo "Running STR genotyping with file set: $fset"

# Set configuration variables; Run STR genotyping 
. ${SCRIPT_DIR}/bpipe_config.sh -t 12 -m 24 -l str_genotype -v $bpipe_ver

MAX_JAVA_MEM=2g ${BPIPE_EXE} run -n $MAX_THREADS -rf "${LOG_REPORT}.html" -p pool="$pool" -p rlen="$rdlen" \
                    -p snvs="$snvs" -p clean="$clean" ${PIPES_DIR}/msi_str_genotyping.pipe ${fset}
