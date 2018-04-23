#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

$0 --- reheader based on provided SQ lines

Version: 20150612

Descriptions:
	reheader BAM
	Index new BAM

Options:
  -h    Print this help message

Example:
  $0 BAM_in SQ_header BAM_out

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
#################### Defaults #######################################



#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
	if command -v $1 >/dev/null 2>&1; then
		exit 0
	else
		echo "I require $1 but it's not installed.  Aborting..." >&2 
		exit 1
	fi
}



###check cmd
#test_samtools=$(CmdExists samtools);
#if (( $(CmdExists samtools) == "1" )); then
#	exit 1
#fi
#test_bamverify=$(CmdExists bamverify);
#if [ $test_bamverify -eq 1 ]; then
#	exit 1
#fi


#################### Defaults #######################################
echo "BAM: $1" >&2
echo "SQ:  $2" >&2
echo "OUT: $3" >&2
bamin=$1
if [ ! -s $bamin ]; then
	echo "Error: input bam  not found" >&2
	exit 1
fi
bamname=${bamin##*/}
bambase=${bamname%.*}
if [ -z "$2" ]; then
	echo "Error: header required" >&2
	exit 1
fi
header=$2
if [ ! -s $header ]; then
	echo "Error: header not found" >&2
	exit 1
fi

if [ -z "$3" ]; then
	bamout="$bambase.reheader.bam"
else
	bamout=$3
fi
if [ -s $bamout ]; then
	echo "Warnings: BAM output exists" >&2
	exit 1
fi
mergefiles=''
echo "###### Summary #####"
echo "Input: $bamin"
echo "Header: $header"
echo "Output: $bamout"
echo -e "\n\n\n"



#################### Input and Output ###############################






#################### Main ###########################################
tempfiles=''
#samtools view -H $bamin > header_ori
#cat header_ori | grep -E "^@HD" > header1
#cat header_ori | grep -E "^@RG" > header2
#cat header_ori | grep -E "^@PG" > header3
#cat header1 $header header2 header3 > header.sam
#samtools reheader header_new $bamin > $bamout
#if [ -s $bamout ]; then
#	rm header_ori header1 header2 header3
#	exit 0
#else
#	echo >&2 "BAM reheader failed: $bamin"
#	exit 1
#fi
samtools view -H $bamin > $bambase.header1
if [ $? -ne 0 ] || [ ! -s $bambase.header1 ]; then
	echo "Error: extract original header" >&2
	exit 1
fi
cat $bambase.header1 | grep -E "^@HD" > $bambase.header2.hd
if [ -s $bambase.header2.hd ]; then
  	mergefiles="$mergefiles $bambase.header2.hd"
  	tempfiles="$tempfiles $bambase.header2.hd"
else
	echo "Error: HD" >&2
	exit 1
fi
mergefiles="$mergefiles $header"
cat $bambase.header1 | grep -E "^@RG" > $bambase.header2.rg
if [ -s $bambase.header2.rg ]; then
  	mergefiles="$mergefiles $bambase.header2.rg"
  	tempfiles="$tempfiles $bambase.header2.rg"
else
	echo "Warnings: RG" >&2
fi
cat $bambase.header1 | grep -E "^@PG" > $bambase.header2.pg
if [ -s $bambase.header2.pg ]; then
	mergefiles="$mergefiles $bambase.header2.pg"
	tempfiles="$tempfiles $bambase.header2.pg $bambase.header1"
else
	echo "Warnings: PG" >&2
fi
cat $mergefiles > $bambase.header2
if [ -s $bambase.header2 ]; then
	samtools view $bamin | cat $bambase.header2 - | samtools view -bSh - > $bamout
	if [ $? -eq 0 ]; then
		bamverify $bamout
		if [ $? -eq 0 ]; then
			if [ ! -s $bamout.bai ]; then
				samtools index 	$bamout
				if [ $? -eq 0 ] && [ -s $bamout.bai ]; then
					rm $tempfiles $bambase.header2
					exit 0
				else
					echo "Index failed" >&2
					exit 1
				fi
			fi
		else
			echo "Nota valid BAM: $bamout" >&2
			exit 1
		fi
	else
		echo "SAMtools reheader failed: $bamin"
		exit 1
	fi
else
	echo "New Header not found: $bamin" >&2
	exit 1
fi

exit 0
