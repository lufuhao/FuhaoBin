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

Version: v20210101

Requirements:
	sed, grep

Descriptions:
	This script is to transform eggnog GO.annotations to input of clusterProfiles

Options:
  -h    Print this help message
  -i    Input eggnog GO.annotations file
  -o    Output

Example:
  $0 -i Split1.annotations -i Split2.annotations -o clusterProfiles.out

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
declare -a AnnotArr=()
opt_o=""
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
#    -i) AnnotArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -o) opt_o=$2;shift 2;;
    -i) AnnotArr=(${AnnotArr[@]} "$2");shift 2;;
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
    echo "Error: CMD $1 not found" >&2
    exit 100
  fi
}





#################### Command test ###################################
CmdExists 'sed'
CmdExists 'grep'



#################### Defaults #######################################




#################### Input and Output ###############################
#emapper.py --cpu 20 -i filter.pep.fa_16 --output filter.pep.fa_16.out -d virNOG  -m diamond
#or
#emapper.py -m diamond \
#           -i sesame.fa \
#           -o diamond \
#           --cpu 19
if [ -e $opt_o ]; then
	rm $opt_o
fi



#################### Main ###########################################
##
#sed -i '/^# /d' diamond.emapper.annotations 
#sed -i 's/#//' diamond.emapper.annotations
#grep ^'#query_name' ${AnnotArr[0]} | sed ' s/^#//' > $opt_o
#numhead_in=$(cat $opt_o | wc -l) 
#if [ $numhead_in -ne 1 ]; then
#	echo "Error: invalid header detected" >&2
#	exit 100
#fi

echo -e "query_name\tseed_eggNOG_ortholog\tseed_ortholog_evalue\tseed_ortholog_score\tPredicted_taxonomic_group\tPredicted_protein_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG Reaction\ttax_scope\teggNOG_OGs\tbestOG\tCOG_Functional_Category\teggNOG_desc" > $opt_o
numhead_in=$(cat $opt_o | wc -l) 




for indAnnot in ${AnnotArr[@]}; do
	numanot=$(grep -v ^'#' $indAnnot | wc -l)
	if [ $numanot -gt 0 ]; then
		echo "File: $numanot Lines: $indAnnot"
		grep -v ^'#' $indAnnot >> $opt_o
		((numhead_in=$numhead_in+$numanot))
	else
		echo "Warnings: File 0 lines: $indAnnot" >&2
	fi
done
echo "Total input line: $numhead_in"
numhead_out=$(cat $opt_o | wc -l)
echo "Total output line: $numhead_out"
if [ $numhead_in -ne $numhead_out ]; then
	echo "Error: Line error: input != output" >&2
	exit 100
fi

exit 0
