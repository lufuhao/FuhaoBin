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

Version: 20180910

Requirements:
	perl
	bedtools

Descriptions:
    * convert BAMPE to BEDPE
    * Modify BAM file to adjust 9 bp for ATAC-seq:
        - forward strand start +4 and
        - reverse strans start -5

Options:
  -h    Print this help message
  -i    Name-sorted BAM
  -o    BEDPE output

Example:
  $0 -i in.bam -o out.bed

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
#[ -z "$@" ] && help
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_i=""

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
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
if CmdExists 'bedtools'; then
	echo "Info: CMD 'bedtools' detected"
else
	echo "Error: CMD 'bedtools' in PROGRAM 'BEDtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if CmdExists 'perl'; then
	echo "Info: CMD 'perl' detected"
else
	echo "Error: CMD 'perl' in PROGRAM 'perl' is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################




#################### Input and Output ###############################
if [ -z "$opt_i" ]; then
	echo "Error: please specify input BAM: -i" >&2
	exit 100
elif [ ! -s "$opt_i" ]; then
	echo "Error: invalid input BAM: -i" >&2
	exit 100
fi
if [ -z "$opt_o" ]; then
	echo "Error: please specify output BED: -o" >&2
	exit 100
fi




#################### Main ###########################################

if [ ! -s "$opt_o.bamtobed" ]; then
	echo "Step(1/2)Info: comvert BAM to bed using bedtools bamtobed"
	bedtools bamtobed -i $opt_i -bedpe > $opt_o.bamtobed 2> $opt_o.bamtobed.err
	if [ $? -ne 0 ] || [ ! -s "$opt_o.bamtobed" ]; then
		echo "Step(1/2)Error: $ProgramName error" >&2
		echo "CMD used: bedtools bamtobed -i $opt_i -bedpe > $opt_o.bamtobed" >&2
		exit 100
	fi
else
	echo "Step(1/2)Info: using existing $opt_o.bamtobed"
fi

if [ ! -s "$opt_o" ]; then
	perl -lane 'if ($F[0] ne $F[3]) {print STDERR "Warnings: read pairs mapped to different chroms: $_"; next;} if ($F[8] eq "+" and $F[9] eq "-") {$F[1]+=4;$F[5]-=5; print "$F[0]\t$F[1]\t$F[5]";}elsif($F[8] eq "-" and $F[9] eq "+"){$F[4]+=4;$F[2]-=5; print "$F[0]\t$F[4]\t$F[2]";}else {print STDERR "Warnings: invalid line: $_";}' $opt_o.bamtobed | sort -k1,1 -k2,2n -k3,3n > $opt_o
	if [ $? -ne 0 ] || [ ! -s $opt_o ]; then
		echo "Step(2/2)Error: $ProgramName error" >&2
		exit 100
	fi
else
	echo "Step(2/2)Info: using existing $opt_o.bamtobed"
fi

exit 0
