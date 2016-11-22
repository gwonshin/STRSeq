#!/bin/bash
# Script sets general parameters needed for bpipe runs 
# Example: bipe_config.sh -t 36 -m 48 -l strseq  (runs with threads=36, mem=48GB, log file prefix='strseq')

TOOLS_DIR=/path/to/tools(bpipe)
bpipe_ver='0.9.8.7'
STRSEQ_DIR='/path/to/strseq'
PIPES_DIR=${STRSEQ_DIR}/Scripts/bpipes

# Turn on extended globbing (allowing regex in file matching patterns)
shopt -s extglob

# Function to print help
#
print_usage()
{
	   echo "Usage: $0 -t {threads} -m {memory(GB)} -l {log prefix} [-v {bpipe version}]"; 
	   echo "Where -t maximum number of threads to use in parallel processing";
	   echo "      -m maximum memory(GB) to use in parallel processing";
	   echo "      -l prefix for bpipe html log file (will be appended with date/time stamp)";
	   echo "      -v bpipe version number";
	   return
}

# Set default parameters   
threads=48
memory=36
logfile=BpipeLog

# Parse command line options
OPTIND=1
while getopts "t:m:l:v:" OPT
do
  case "$OPT" in
    t) threads=$OPTARG;;
    m) memory=$OPTARG;;
    l) logfile="$OPTARG";;
	v) bpipe_ver="$OPTARG";;
   \?) print_usage; exit 1;;
  esac
done
	
if [ $threads -lt 1 ] || [ $threads -gt 60 ]; then
  echo "Error: Threads must be between 1 and 60"; exit 2
fi

if [ $memory -lt 4 ] || [ $memory -gt 156 ]; then
  echo "Error: Memory must be between 4 and 156"; exit 3
fi

BPIPE_EXE=${TOOLS_DIR}/bpipe-${bpipe_ver}/bin/bpipe
if [ ! -f $BPIPE_EXE ]; then
  echo "Error: bpipe executable not found at: ${BPIPE_EXE}"; exit 4
fi

MAX_THREADS=$threads
GB_MEM=$memory
MB_MEM=$((memory*1024))

timestamp=$( date +%y%m%d%H%M)
LOG_REPORT="${logfile}_${timestamp}"

echo "Threads: $MAX_THREADS, Memory(MB): $MB_MEM, Log file: $LOG_REPORT"
echo "Bpipe: ${BPIPE_EXE}"

