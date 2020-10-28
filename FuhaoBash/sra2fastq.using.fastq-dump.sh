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

$0 --- Convert SRA to fastq.gz

Version: v20200723

Requirements:
	SRAtoolkit

Descriptions:
	Convert SRA to fastq files

Options: fastq-dump options
  -h    Print this help message
  -i    SRA file, comma delimited
  -d    Ouput Path
  -n    Fastq-dump options
          Defaults: --gzip --split-3

Example:
  $0 -i SRR000000.sra,SRR111111.sra

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
opt_d=$PWD
declare -a SraArr=()
opt_n=' --gzip --split-e '
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) SraArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -d) opt_d=$2;shift 2;;
    -n) opt_n=$2;shift 2;;
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
CmdExists 'fastq-dump'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'fastq-dump' in PROGRAM 'SRAtoolkit' is required but not found.  Aborting..." >&2 
	exit 127
fi

#################### Defaults #######################################




#################### Input and Output ###############################




#################### Main ###########################################

for ind_sra in ${SraArr[@]}; do
	fastq-dump $opt_n --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' --outdir $opt_d $ind_sra
	if [ $? -ne 0 ]; then
		echo "Error: fastq-dump running error" >&2
		echo "CMD used: fastq-dump $opt_n --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' --outdir $opt_d $ind_sra" >&2
		exit 100;
	fi
done

exit 0
