#!/bin/bash
### Source: http://data.biostarhandbook.com/scripts/wonderdump.sh
set -ue
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

$0 --- download SRR using curl

Version: v20200723

Requirements:
	curl

Descriptions:
  A workaround to download SRA files directly using curl
    when fastq-dump's internet connection does not work.
    Which can happen surprisingly frequently.

Options: 
  -h    Print this help message
  -f    SRR list file
  -i    SRR names, comma delimited
  -d    Output Path

Example:

  $0 -f SRR_list.file -i SRRxxxxxx,SRRyyyyyyy -d $HOME/puclic/SRA

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
declare -a SRRlist=()
opt_f=""
opt_d=$PWD

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) SRRlist=($(echo $2 | tr "," "\n"));shift 2;;
    -f) opt_f=$(readlink -f $2);shift 2;;
    -d) opt_d=$2;shift 2;;
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
SrrDownloadUsingCurl() {
	local $SDUCsrr=$1
	
	local SDUCsub="(SrrDownloadUsingCurl)"
	# Create the full path to the file.
	local Sra_File="$opt_d/$SDUCsrr.sra"
	local Tmp_File="$opt_d/$SDUCsrr.tmp"
	
	echo "($SDUCsub)Info: Getting SRR run: $SDUCsrr"

	cd $opt_d
	# Download only if it does not exist.
	if [ ! -f $Sra_File ]; then
		PATH1=${SDUCsrr:0:6}
		PATH2=${SDUCsrr:0:10}
		URL="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${PATH1}/${PATH2}/${SDUCsrr}.sra"
		echo "($SDUCsub)Info: Downloading: $URL"
		echo "($SDUCsub)Info: Saving to: $Sra_File"
		curl $URL > $Tmp_File 2> $Sra_File.curl.err
		if [ $? -ne 0 ]; then
			echo "($SDUCsub)Error: curl failed to download SRR: $SDUCsrr" >&2
			return 100;
		else # Move to local file only if successful.
			mv $Tmp_File $Sra_File
			if [ -e $Sra_File.curl.err ]; then
				rm $Sra_File.curl.err >/dev/null 2>&1
			fi
		fi
	else
		echo "($SDUCsub)Warnings: existing SRA file found: $Sra_File"
	fi

	return 0
	
}



#################### Command test ###################################
CmdExists 'curl'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'curl' in PROGRAM 'curl' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'fastq-dump'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'fastq-dump' in PROGRAM 'SRAtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################



#################### Input and Output ###############################
# Make the directory if it does not exist.
if [ ! -d $opt_d ]; then
	mkdir -p $opt_d
fi



#################### Main ###########################################

cd $opt_d
if [ ! -z "$opt_f" ] && [ -s "$opt_f" ]; then
	echo ""
	echo "### Using SRR list file: $opt_f"
	echo ""
	while read SrrID; do
		echo "###     SRR: $SrrID"
		if SrrDownloadUsingCurl $SrrID; then
			echo "###     SRR: $SrrID    OK"
		else
			echo "###     SRR: $SrrID    failed" >&2
		fi
	done < $opt_f
fi

if [[ ${SRRlist[@]} -gt 0 ]]; then
	echo ""
	echo "### Using SRR comma list"
	echo ""
	for SrrID in ${SRRlist[@]}; do
		echo "###     SRR: $SrrID"
		if SrrDownloadUsingCurl $SrrID; then
			echo "###     SRR: $SrrID    OK"
		else
			echo "###     SRR: $SrrID    failed" >&2
		fi
	done
fi

exit 0
