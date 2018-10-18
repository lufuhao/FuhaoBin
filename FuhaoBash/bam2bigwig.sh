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

Version: 20181010

Requirements:
  bedGraphToBigWig
    http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
  bedtools
  samtools
  bam_clean_header.pl

Descriptions:
    convert BAM to bedGraph and then bigwig
      1. extract BAM if -r specified
           samtools
           bam_clean_header.pl
      2. bam to bedgraph if -i specified
           bedtools
           Will ignore this step if specified -bg but not -i
      3. bedGraph to bigwig
           bedGraphToBigWig

Options:
  -h    Print this help message
  -i    BAM file
  -fai  genome fasta index input
  -bg   bedGraph output [BAMbase.bedGraph] or input [existing bedgraph file]
  -r    Special regions: Chr1:1-1000
  -o    bigwig output

Example:
  $0 -i my.bam -fai myGenome.fa.fai -o 

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
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n";
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_i="";
opt_fai="";
opt_bg="";
opt_r="";
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -r) opt_r=$2;shift 2;;
    -fai) opt_fai=$2;shift 2;;
    -bg) opt_bg=$2;shift 2;;
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
CmdExists 'bedGraphToBigWig'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bedGraphToBigWig' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bedtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bedtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ ! -z "$opt_r" ]; then
	CmdExists 'samtools'
	if [ $? -ne 0 ]; then
		echo "Error: CMD 'samtools' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	CmdExists 'bam_clean_header.pl'
	if [ $? -ne 0 ]; then
		echo "Error: CMD 'bam_clean_header.pl' is required but not found.  Aborting..." >&2 
		exit 127
	fi
fi


#################### Defaults #######################################
RED='\033[0;31m'
NC='\033[0m' # No Color
#printf "I ${RED}love${NC} Stack Overflow\n"


#################### Input and Output ###############################

if [ -z "$opt_fai" ]; then
	echo "Error: need -fai" >&2
	exit 100
elif [ ! -s $opt_fai ]; then
	echo "Error: invalid fai" >&2
	exit 100
fi




#################### Main ###########################################
#you can create a bedgraph directly from bam with bedtools alone:

if [ -z "$opt_i" ]; then
	echo "Step1: Region not specified, ignore..."
	if [ -z "$opt_bg" ]; then
		echo "Error: invalid input; specify either -i or -bg" >&2
		exit 100
	elif [ -s $opt_bg ]; then
		echo -e "${RED}Warnings: ignore BAM to bedGraph${NC}" >&2
		echo "Step2: using existsing bedGFraph file as input: $opt_bg"
		if [ -z "$opt_o" ]; then
			opt_o="$opt_bg.bw"
		fi
	fi
elif [ -s $opt_i ]; then
	MyBamName=${opt_i##*/}
	if [ ! -z "$opt_r" ]; then
		echo "Step1: Extracting Region: $opt_r"
		BamFilter="$MyBamName.filter.bam"
#		echo "samtools view -b -h $opt_i $opt_r > $BamFilter"
		samtools view -b -h $opt_i $opt_r > $BamFilter
		if [ $? -ne 0 ] || [ ! -s $BamFilter ]; then
			echo "Step1 error: bedtools genomecov" >&2
			echo "CMD used: samtools view -b -h $opt_i $opt_r > $BamFilter" >&2
			exit 100
		fi
		BamFilter2="$MyBamName.filter2.bam"
		bam_clean_header.pl -i $BamFilter -o $BamFilter2
		if [ $? -ne 0 ] || [ ! -s $BamFilter2 ]; then
			echo "Step1 error: bam_clean_header.pl" >&2
			echo "CMD used: bam_clean_header.pl -i $BamFilter -o $BamFilter2" >&2
			exit 100
		fi
		rm $BamFilter > /dev/null 2>&1
		opt_i=$BamFilter2
	fi
	if [ -z "$opt_bg" ]; then
		opt_bg="$MyBamName.bedGraph"
	fi
	if [ -z "$opt_o" ]; then
		opt_o="$opt_bg.bw"
	fi
	echo "Step2: convert BAM to bedGraph"
	echo "       BAM input: $opt_i"
	echo "       BedGraph out: $opt_bg"
	bedtools genomecov -bga -ibam $opt_i > $opt_bg
	if [ $? -ne 0 ] || [ ! -s $opt_bg ]; then
		echo "Step2 error: bedtools genomecov" >&2
		echo "CMD used: bedtools genomecov -bga -ibam $opt_i > $opt_bg" >&2
		exit 100
	fi
else
	echo -e "${RED}Warnings: invalid BAM input${NC}" >&2
	exit 100
fi


#you can then convert it to bigwig directly:
echo "Step3: convert bedGraph to bigWig"
bedGraphToBigWig $opt_bg $opt_fai $opt_o
if [ $? -ne 0 ] || [ ! -s $opt_o ]; then
	echo "${RED}Error: Step3 bedGraphToBigWig${NC}" >&2
	echo "CMD used: bedGraphToBigWig $opt_bg $opt_fai $opt_o" >&2
	exit 100
fi
