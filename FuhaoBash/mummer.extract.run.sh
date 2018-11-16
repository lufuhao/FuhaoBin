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

Version: 20171219

Requirements:
	Script: mum.stat
	CMD: seqtk

Descriptions:
	Extract seqids from multifasta
	Run mum.stat

Options:
  -h    Print this help message
  -r    Reference database
  -q    Query database
  -1    Reference seqids: accept multi-value
  -2    Query seqids: accept multi-value
  -p    Output prefix

Example:
  $0 -r ref.fa -q qry.fa -1 rfn1 -1 rfn2 -2 qry1 -2 qry2 -p MyTemp

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
declare -a seq_rfn
declare -a seq_qry
opt_r=''
opt_q=''
opt_p='MyTemp'
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -r) opt_r=$2;shift 2;;
    -q) opt_q=$2;shift 2;;
    -1) seq_rfn=(${seq_rfn[@]} "$2");shift 2;;
    -2) seq_qry=(${seq_qry[@]} "$2");shift 2;;
    -p) opt_p=$2;shift 2;;
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
### Extract fasta by a file of read ID list
SeqTkSubSeqFasta () {
	local SSFfastain=$1;
	local SSFidlist=$2;
	local SSFfastaout=$3;
	local SSFsubinfo='SUB(SeqTkSubSeqFasta)';
	
	local SSFnumberlist=0;
	local SSFnumberout=0;
	local SSFfa_line_width=70;
	
	if [ -z "$SSFfastain" ] || [ ! -s "$SSFfastain" ]; then
		echo "${SSFsubinfo}Error: invalid input fasta: $SSFfastain" >&2
		return 1;
	fi
	if [ -z "$SSFidlist" ] || [ ! -s "$SSFidlist" ]; then
		echo "${SSFsubinfo}Error: invalid fasta ID list: $SSFidlist" >&2
		return 1;
	fi
	if [ -z "$SSFfastaout" ]; then
		echo "${SSFsubinfo}Error: invalid fasta ID list: $SSFfastaout" >&2
		return 1;
	fi
	if [ -e $SSFfastaout ]; then
		rm -rf $SSFfastaout >/dev/null 2>/dev/null
	fi
	
	seqtk subseq -l $SSFfa_line_width $SSFfastain $SSFidlist > $SSFfastaout
	if [ $? -ne 0 ] || [ ! -s $SSFfastaout ]; then
		echo "${SSFsubinfo}Error: seqtk subseq fasta error" >&2
		echo "${SSFsubinfo}CMD used: seqtk subseq -l 70 $SSFfastain $SSFidlist > $SSFfastaout" >&2
		return 1;
	fi
	
	SSFnumberlist=$(wc -l $SSFidlist)
	SSFnumberout=$(perl -ne 'BEGIN{$linnum=0;} $linenum++ if (/^>/); END {print $linenum, "\n";}' < $SSFfastaout)

	if [ -z "$SSFnumberlist" ] || [ $SSFnumberlist -eq 0 ] || [ $SSFnumberlist -ne $SSFnumberout ]; then
		echo "${SSFsubinfo}Error: seqtk lineout partially failed" >&2
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist" >&2
		echo "${SSFsubinfo}        Total number of R1 IDs extracted: $SSFnumberout" >&2
		return 1;
	else
		echo "${SSFsubinfo}Info: seqtk subseq secceeded"
		echo "${SSFsubinfo}        Total number of IDs to extract:   $SSFnumberlist"
	fi

	return 0;
}




#################### Command test ###################################
if [[ $(CmdExists 'mum.stat') -eq 1 ]]; then
	echo "Error: script 'mum.stat' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [[ $(CmdExists 'seqtk') -eq 1 ]]; then
	echo "Error: CMD 'seqtk' in PROGRAM 'seqtk' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
mum_identity=95
list_rfn_ids="$opt_p.rfn"
list_rfn_fasta="$list_rfn_ids.fa"
list_qry_ids="$opt_p.qry"
list_qry_fasta="$list_rfn_ids.fa"



#################### Input and Output ###############################
if [ -z "$opt_r" ] || [ ! -s "$opt_r" ]; then
	echo "Error: invalid reference\n" >&2
	exit 100
fi
if [ -z "$opt_q" ] || [ ! -s "$opt_q" ]; then
	echo "Error: invalid query\n" >&2
	exit 100
fi
if [[ ${#seq_rfn[@]} -eq 0 ]]; then
	echo "Error: invalid reference list\n" >&2
	exit 100
fi
if [[ ${#seq_qry[@]} -eq 0 ]]; then
	echo "Error: invalid query list\n" >&2
	exit 100
fi


#################### Main ###########################################

echo "Reference:  $opt_r"
echo "    list:   $list_rfn_ids"
echo "    fasta:  $list_rfn_fasta"
echo "Query:      $opt_q"
echo "    list:   $list_qry_ids"
echo "    fasta:  $list_qry_fasta"
echo "Total ${#seq_rfn[@]} reference seqids"
echo "Total ${#seq_qry[@]} query seqids"

num_qry=0;
num_rfn=0;
if [ -e "$list_rfn_ids" ] || [ -e "$list_rfn_fasta" ]; then
	rm "$list_rfn_ids" "$list_rfn_fasta" > /dev/null 2>&1
fi
if [ -e "$list_qry_ids" ] || [ -e "$list_qry_fasta" ]; then
	rm "$list_qry_ids" "$list_qry_fasta" > /dev/null 2>&1
fi
###
echo "Extract reference IDs"
for ind_rfn in "${seq_rfn[@]}"; do
	echo "$ind_rfn"
	echo "$ind_rfn" >> $list_rfn_ids
	((num_rfn++));
	echo "Total $num_rfn reference seqids"
done
echo "Total $num_rfn reference seqids"
if [ $num_rfn -gt 0 ]; then
	echo "Total $num_rfn reference seqids"
	for ind_rfn in ${seq_rfn[@]}; do
		echo "        $ind_rfn"
	done
else
	echo "Error: no reference IDs detected" >&2
	exit 100
fi
echo "Extract reference seqs"
if SeqTkSubSeqFasta $opt_r $list_rfn_ids $list_rfn_fasta; then
	echo "Info: Reference ids extracted"
else
	echo "Error: failed to extract reference seqs" >&2
	exit 100
fi
###
echo "Extract query IDs"
for ind_qry in ${seq_qry[@]}; do
	echo "$ind_qry" >> $list_qry_ids
	((num_qry++));
done
echo "Total $num_qry query seqids"
if [ $num_qry -gt 0 ]; then
	echo "Total $num_qry query seqids"
	for ind_qry in ${seq_qry[@]}; do
	echo "        $ind_qry"
	done
else
	echo "Error: no query IDs detected" >&2
	exit 100
fi
if SeqTkSubSeqFasta $opt_q $list_qry_ids $list_qry_fasta; then
	echo "Info: query ids extracted"
else
	echo "Error: failed to extract query seqs" >&2
	exit 100
fi

###
if [ -s $list_rfn_fasta ] && [ -s $list_qry_fasta ]; then
	echo "Running mummer"
	mum.stat -i $mum_identity -r $list_rfn_fasta -q $list_qry_fasta -o $opt_p
	if [ $? -ne 0 ]; then
		echo "Error: mum.stat running" >&2
		exit 100
	fi
	if [ -s "$opt_p.filter.q" ]; then
		show-coords -c -r -l -T $TEMPID.temp.filter.q
	fi
else
	echo "Error: mum.stat input files" >&2
	exit 100
fi


exit 0;
