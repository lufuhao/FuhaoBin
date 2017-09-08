#!/bin/bash
### Exit if command fails
set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
#set -o nounset
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

#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
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

Version: 20170904

Requirements:
	perl && File::Spec

Descriptions:
	Merge nucmer-generated delta files

Options:
  -h    Print this help message
  -i    Comma-delimited delta files
  -f    File containing with list, one file per line
  -r    Full path to Reference fasta
  -q    Full path to Query fasta
  -o    File for final output delta
  -debug

Example:
  $0 -i part1.delta,part2delta -r $HOME/ref.fa -q $HOME/qry.fa

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
opt_i='';
opt_f='';
opt_r='';
opt_q='';
opt_o='MyOutput.final.delta';
opt_debug=0;
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -f) opt_f=$2;shift 2;;
    -r) opt_r=$2;shift 2;;
    -q) opt_q=$2;shift 2;;
    -o) opt_o=$2;shift 2;;
    -debug) opt_debug=1;shift 1;;
    --) shift;break;;
    -*) echo "Error: no such option $1. -h for help" > /dev/stderr;exit 1;;
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
if [[ $(CmdExists 'perl') -ne 0 ]]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
declare -a file1list=()
declare -a file2list=()
declare -a file3list=()


#################### Input and Output ###############################
if [[ "$opt_r" =~ \/ ]] && [ -s "$opt_r" ] && [[ "$opt_q" =~ \/ ]] && [ -s "$opt_q" ]; then
	echo "Info: Ref: $opt_r"
	echo "Info: Qry: $opt_q"
else
	opt_r=$(echo $(cd $(dirname "$opt_r"); pwd)/$(basename "$opt_r"))
	opt_q=$(echo $(cd $(dirname "$opt_q"); pwd)/$(basename "$opt_q"))
	if [[ "$opt_r" =~ \/ ]] && [ -s "$opt_r" ] && [[ "$opt_q" =~ \/ ]] && [ -s "$opt_q" ]; then
		echo "Info: Ref: $opt_r"
		echo "Info: Qry: $opt_q"
	else
		echo "Error: please use full path of reference and query" >&2
		exit 100;
	fi
fi



#################### Main ###########################################



if [ "$opt_i" == '' ]; then
	echo "Info: Comma-delemited list -i NOT detected"
else
	echo "Info: Comma-delemited list -i detected"
	file1list=($(echo "$opt_i" | perl -lne '@arr=split(/,/); foreach $x (@arr) {print $x;}'))
fi
if [ "$opt_f" == '' ]; then
	echo "Info: delta list file -f NOT detected"
else
	echo "Info: delta list file -f detected"
	file2list=($(perl -lne 'print;' < $opt_f))
fi

file3list=("${file1list[@]}" "${file2list[@]}")
echo "Info: total number of delta: "${#file3list[@]}
echo "Info: total number of delta: "${#file3list[@]} >&2

for mydelta in "${file3list[@]}"; do
	if [ -s "$mydelta" ]; then
		if [[ $opt_debug -eq 1 ]]; then
			echo "    existing    $mydelta"
		fi
	else
		echo "    nonexist    $mydelta" >&2
		exit 100
	fi
done

echo "$opt_r $opt_q" > $opt_o
echo "NUCMER" >> $opt_o

for mydelta in "${file3list[@]}"; do
	Echo 
	tail -n +3 $mydelta >> $opt_o
done

if [ -s $opt_o ]; then
	exit 0;
else
	echo "Error: Delta merge Failed" >&2
	exit 100
fi
