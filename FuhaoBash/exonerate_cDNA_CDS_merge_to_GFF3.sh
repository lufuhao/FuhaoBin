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

$0 --- merge Exonerate-mapped cDNA and CDS into GFF3

Version: 20180227

Requirements:
	Linux: grep, perl, cat
	Script: exonerate_GFF3_reformat.pl

Descriptions:
	[Separately Map cDNA and CDS onto references using Exonerate]
	Merge cDNA and CDS into GFF3
		* Will generate temporary files: -1 -2 -3; finally delete

Options:
  -h    Print this help message
  -i    GFF3 gene list file
  -n    cDNA GFF3
  -s    CDS GFF3
  -o    output GFF3
  -1    Temporary cDNA file [TEMP__CDNA.gff3]
  -2    Temporary CDS file [TEMP__CDS.gff3]
  -3    Temporary merge GFF3 file [TEMP__MERGE.gff3]

Example:
  $0 -i my.gene.list -n cDNA.gff3 -s CDS.gff3 -o output.gff3

  Gene list
  grep -P "\tgene\t" xxx.gff3 | cut -f 9 | perl -lane 's/^.*ID=//; s/;.*$//; print;' | sort -u > genelist

  * test if there are repeats
  grep -P "\tmRNA\t" wheat2tauschii.list.gff3 | cut -f 9 | perl -lne 's/^.*ID=//;s/;.*$//;print;' | sort | uniq -cd

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
opt_i=''
opt_n=''
opt_s=''
tempfile_cdna='TEMP__CDNA.gff3'
tempfile_cds='TEMP__CDS.gff3'
tempfile_merge='TEMP__MERGE.gff3'
debug=0
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -n) opt_n=$2;shift 2;;
    -s) opt_s=$2;shift 2;;
    -o) opt_o=$2;shift 2;;
    -1) tempfile_cdna=$2;shift 2;;
    -2) tempfile_cds=$2;shift 2;;
    -3) tempfile_merge=$2;shift 2;;
    -debug) debug=1; shift 1;;
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
CleanTempFiles () {
	if [ $debug -ne 0 ]; then
		echo "    * Cleaning temporary files"
	fi
	if [ -e "$tempfile_cdna" ]; then
		rm "$tempfile_cdna" > /dev/null 2>&1
	fi
	if [ -e "$tempfile_cds" ]; then
		rm "$tempfile_cds" > /dev/null 2>&1
	fi
	if [ -e "$tempfile_merge" ]; then
		rm "$tempfile_merge" > /dev/null 2>&1
	fi
	return 0
}

#################### Command test ###################################
if [[ $( CmdExists 'exonerate_GFF3_reformat.pl' ) -ne 0 ]]; then
	echo "Error: script 'exonerate_GFF3_reformat.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################
if [ -z "$opt_i" ] || [ ! -s "$opt_i" ]; then
	echo "Error: invalid Gene List file for option '-i'" >&2
	exit 100
fi
if [ -z "$opt_n" ] || [ ! -s "$opt_n" ]; then
	echo "Error: invalid cDNA GFF3 file for option '-n'" >&2
	exit 100
fi
if [ -z "$opt_s" ] || [ ! -s "$opt_s" ]; then
	echo "Error: invalid CDS GFF3 file for option '-s'" >&2
	exit 100
fi
if [ -e "$opt_o" ]; then
	echo "Error: existing output GFF3 file for option '-o'" >&2
	exit 100
fi
if [ -z "$tempfile_cdna" ]; then
	echo "Error: invalid cDNA GFF3 file for option '-1'" >&2
	exit 100
fi
if [ -z "$tempfile_cds" ]; then
	echo "Error: invalid CDS GFF3 file for option '-2'" >&2
	exit 100
fi
if [ -z "$tempfile_merge" ]; then
	echo "Error: invalid merge GFF3 file for option '-3'" >&2
	exit 100
fi




#################### Main ###########################################
CleanTempFiles

while read GeneName; do
	
	echo "##gff-version 3" > "$tempfile_cdna"
	echo "##gff-version 3" > "$tempfile_cds"
	
	echo "Gene: $GeneName"
	
	declare -a genenames1=()
	genenames1=($(grep "$GeneName" "$opt_n" | perl -lne 's/^.*ID=//;s/;.*$//;print;'))
	if [ ${#genenames1[@]} -gt 0 ]; then
		for GeneID1 in ${genenames1[@]}; do
			echo "    Transferable    cDNA    $GeneID1"
			grep "$GeneID1" "$opt_n" >> "$tempfile_cdna"
			
		done
	else
		echo "    Transferable    cDNA   NONE" >&2
	fi

	declare -a genenames2=()
	genenames2=($(grep "$GeneName" "$opt_s" | perl -lne 's/^.*ID=//;s/;.*$//;print;'))
	if [ ${#genenames2[@]} -gt 0 ]; then
		for GeneID2 in ${genenames2[@]}; do
			echo "    Transferable    CDS    $GeneID2"
			grep "$GeneID2" "$opt_s" >> "$tempfile_cds"
		done
	else
		echo "    Transferable    CDS   NONE" >&2
	fi
	
	exonerate_GFF3_reformat.pl "$tempfile_cdna" "$tempfile_cds" "$tempfile_merge"
	if [ $? -ne 0 ] || [ ! -s "$tempfile_merge" ]; then
		echo "Error: exonerate_GFF3_reformat.pl running error" >&2
		exit 1
	else
		cat "$tempfile_merge" >> "$opt_o"
	fi
	
	CleanTempFiles

done < "$opt_i"


exit 0
