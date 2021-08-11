#!/bin/bash
### Exit if command fails
#set -o errexit
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

###echo color
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
#RED='\033[0;31m'
#NC='\033[0m' # No Color
#printf "I ${RED}love${NC} Stack Overflow\n"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: v20210811

Requirements:
	perl

Descriptions:
	Given a list of BAMs, use samtools flagstat to do statistics

Options:
  -h    Print this help message
  -i    Individual BAM files, -i bam1 -i bam2
  -f    BAM list file
  -o    Statistics Output, default: "samtools.flagstat.tab"
  -t    Number of threads, default: 1

Example:
  $0 -i ./chr1.fa -i 1.bam -i 2.bam -f bam.list -o samtools.flagstat.tab

Author:
  Fu-Hao Lu
  Professor, PhD
  State Key Labortory of Crop Stress Adaptation and Improvement
  College of Life Science
  Jinming Campus, Henan University
  Kaifeng 475004, P.R.China
  E-mail: lufuhao@henu.edu.cn
HELP
exit 2
}
[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
declare -a bamlist1=()
declare -a bamlist2=()
declare -a bamAll=()
opt_t=1
opt_o=""
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) bamlist1=(${bamlist1[@]} "$2");shift 2;;
    -f) opt_f=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -o) opt_o=$2;shift 2;;
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
CmdExists 'samtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
sam_flagstat_option=""



#################### Input and Output ###############################
if [ ! -z "$opt_f" ]; then
	if [ -s "$opt_f" ]; then
		mapfile bamlist2 < "$opt_f"
	fi
fi
bamAll=("${bamlist1[@]}" "${bamlist2[@]}")


if [ -e "$opt_o.log" ]; then
	echo "Error: existing output: $opt_o.log" >&2
	exit 100
fi
if [ -e "$opt_o" ]; then
	echo "Error: existing output: $opt_o" >&2
	exit 100
fi

if [[ "$opt_t" =~ [0-9]+ ]]; then
	sam_flagstat_option=" --threads $opt_t "
fi

#################### Main ###########################################

for indbam in "${bamAll[@]}"; do
	if [ ! -e "$indbam" ]; then
		echo "Error: BAM not found: $indbam" >&2
		continue
	fi
	echo "$indbam"
	echo "$indbam" >> "$opt_o.log"
	echo "$indbam" >&2
	samtools flagstat $sam_flagstat_option "$indbam" >> "$opt_o.log"
	if [ $? -ne 0 ]; then
		echo "Error: BAM flagstat failed: $indbam" >&2
		continue
	fi
done

perl -ne 's/\s*in total\s*\s*\(QC-passed\s*reads\s*\+\s*QC-failed\s*reads\)\n/\t/; s/\s*secondary\n/\t/;s/\s*supplementary\n/\t/;s/\s*duplicates\n/\t/;s/\s*mapped\s*\(/\t/;s/\s*:\s*N\/A\)\n/\t/;s/\s*paired\s*in\s*sequencing\n/\t/;s/\s*read\d+\n/\t/;s/\s*properly\s*paired\s*\(/\t/;s/\s*with\s*itself\s*and\s*mate\s*mapped\n/\t/;s/\s*singletons\s*\(/\t/;s/\s*with\s*mate\s*mapped\s*to\s*a\s*different\s*chr\s*\(mapQ>=5\)$//;s/\s*with\s*mate\s*mapped\s*to\s*a\s*different\s*chr\n/\t/; s/\s*\+\s*0//g;print;' "$opt_o.log" > "$opt_o"

exit 0
