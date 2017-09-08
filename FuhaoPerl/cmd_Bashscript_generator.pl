#!/usr/bin/env perl
use warnings;
use strict;
###Version: 20170619
print '#!/bin/bash
### Exit if command fails
set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
set -o nounset
### Script Root
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
### MachType
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

#export NUM_THREADS=`grep -c \'^processor\' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20170619

Requirements:
	perl && File::Spec

Descriptions:
	xxx

Options:
  -h    Print this help message
  -i    CONFIG file
  -t    Number of threads, default: 1
  -s    Not run simulation
  -a    Not run assembly

Example:
  $0 -i ./chr1.fa -t 10

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) opt_s=1;shift 1;;
    -a) opt_a=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    return 1
  fi
}





#################### Command test ###################################
if [ $(CmdExists \'santools\') -eq 0 ]; then
	exit 0
else
	echo "Error: CMD/script \'samtools\' in PROGRAM \'SAMtools\' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists \'santools\') -eq 1 ]; then
	echo "Error: CMD/script \'samtools\' in PROGRAM \'SAMtools\' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################




#################### Main ###########################################
if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	exit 1
fi
';
