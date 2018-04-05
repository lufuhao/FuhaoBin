#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
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

Version: 20161025

Requirements:
	Package: bowtie, samtools
	Script: picardrmdup

Descriptions:
	Map each mate to reference separately using 'bowtie'
	Merge two BAMs using 'samtools merge'
	picard rmdup to remove duplicates

Options:
  -h                     Print this help message
  -1  [1.R1.fq]          First mate
  -2  [1.R2.fq]          Second mate
  -b  [Msg]              bowtie options
  -g  [readgroup]        Read group
  -p  [Opt]              Output prefix, default: MyOutput
  -d                     Delete temporary files

Note:
	Bowtie option
	For fastq: ' -q -p 10 -v 0 -a --sam /path/to/index '
	For Fasta: ' -f -p 10 -v 0 -a --sam /path/to/index '

Temporary Files
	$opt_p.1.st.bam    $opt_p.1.st.bam.bai
	$opt_p.2.st.bam    $opt_p.2.st.bam.bai
	$opt_p.1.2.st.merge.bam    $opt_p.1.2.st.merge.bam.bai
	$opt_p.1.2.st.merge.rmdup.bam    $opt_p.1.2.st.merge.rmdup.bam.bai

Example:
  $0 -1 R1.fq -2 R2.fq -b ' -q -p 10 -v 0 -a --sam /path/to/index ' 

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
opt_1=''
opt_2=''
opt_b=''
opt_p='MyOutput'
path_bowtie='bowtie'
path_bowbuild='bowtie-build'
path_samtools='samtools'
path_picardrmdup='picardrmdup'
path_seqtk='seqtk'
tempdel=0
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) opt_1=$2;shift 2;;
    -2) opt_2=$2;shift 2;;
    -d) tempdel=1; shift 1;;
    -b) opt_b=$2;shift 2;;
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
    echo 0
  else
#    echo "I require $1 but it's not installed.  Aborting." >&2
    echo 1
  fi
#  local cmd=$1
#  if command -v $cmd >/dev/null 2>&1;then
#    echo >&2 $cmd "  :  "`command -v $cmd`
#    exit 0
#  else
#    echo >&2 "Error: require $cmd but it's not installed.  Exiting..."
#    exit 1
#  fi
}

###Usage: array=(`split delimiter string`)
split () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}

#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}
abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }


#################### Command test ###################################
#if [ $(CmdExists "$path_bowbuild") -ne 0 ]; then
#	echo "Error: CMD 'bowtie-build' in PROGRAM 'BOWTIE' is required but not found.  Aborting..." >&2 
#	exit 127
#fi
#echo "INFO: bowtie-build path: $path_bowbuild"
if [ $(CmdExists "$path_bowtie") -ne 0 ]; then
	echo "Error: CMD 'bowtie' in PROGRAM 'BOWTIE' is required but not found.  Aborting..." >&2 
	exit 127
fi
echo "INFO: bowtie path: $path_bowtie"
if [ $(CmdExists "$path_samtools") -ne 0 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
echo "INFO: samtools path: $path_samtools"
if [ $(CmdExists "$path_picardrmdup") -ne 0 ]; then
	echo "Error: script 'picardrmdup' is required but not found.  Aborting..." >&2 
	exit 127
fi
echo "INFO: picardrmdup path: $path_picardrmdup"
#if [ $(CmdExists "$path_seqtk") -ne 0 ]; then
#	echo "Error: CMD 'seqtk' in PROGRAM 'seqtk' is required but not found.  Aborting..." >&2 
#	exit 127
#fi
#echo "INFO: seqtk path: $path_seqtk"

#################### Defaults #######################################



#################### Input and Output ###############################
if [ -z "$opt_1" ] || [ ! -s "$opt_1" ]; then
	echo "Error: invalid first mate file" >&2 
	exit 1
fi
if [ -z "$opt_2" ] || [ ! -s "$opt_2" ]; then
	echo "Error: invalid second mate file" >&2 
	exit 1
fi



#################### Main ###########################################
echo "INFO: first mate: $opt_1"
echo "INFO: second mate: $opt_2"
echo "INFO: bowtie options: $opt_b"
echo "INFO: output prefix: $opt_p"



### Bowtie mapping
echo "###Bowtie mapping $opt_p"
echo "###Bowtie mapping $opt_p" >&2
echo "##Bowtie mapping $opt_p.1"
echo "##Bowtie mapping $opt_p.1" >&2
$path_bowtie $opt_b $opt_1  | $path_samtools view -b -h -S -F 4 - > $opt_p.1.bam
if [ $? -ne 0 ] || [ ! -s $opt_p.1.bam ]; then
	echo "Error: bowtie mapping error: $opt_p.1.bam" >&2
	exit 1
fi
echo "##Bowtie mapping $opt_p.2"
echo "##Bowtie mapping $opt_p.2" >&2
$path_bowtie $opt_b $opt_2  | $path_samtools view -b -h -S -F 4 - > $opt_p.2.bam
if [ $? -ne 0 ] || [ ! -s $opt_p.2.bam ]; then
	echo "Error: bowtie mapping error: $opt_p.2.bam" >&2
	exit 1
fi



### Samtools sort
echo "###Samtools sort $opt_p"
echo "###Samtools sort $opt_p" >&2
echo "##Samtools sort $opt_p.1"
echo "##Samtools sort $opt_p.1" >&2
$path_samtools sort "$opt_p.1.bam" "$opt_p.1.st"
if [ $? -ne 0 ] || [ ! -s "$opt_p.1.st.bam" ]; then
	echo "Error: Samtools sort error: $opt_p.1.st.bam" >&2
	exit 1
else
	rm $opt_p.1.bam
fi
echo "##Samtools sort $opt_p.2"
echo "##Samtools sort $opt_p.2" >&2
$path_samtools sort "$opt_p.2.bam" "$opt_p.2.st"
if [ $? -ne 0 ] || [ ! -s "$opt_p.2.st.bam" ]; then
	echo "Error: Samtools sort error: $opt_p.2.st.bam" >&2
	exit 1
else
	rm $opt_p.2.bam
fi


### Samtools index
echo "###Samtools index $opt_p"
echo "###Samtools index $opt_p" >&2
echo "##Samtools index $opt_p.1"
echo "##Samtools index $opt_p.1" >&2
$path_samtools index "$opt_p.1.st.bam"
if [ $? -ne 0 ] || [ ! -s "$opt_p.1.st.bam.bai" ]; then
	echo "Error: Samtools index error: $opt_p.1.st.bam" >&2
	exit 1
fi
echo "##Samtools index $opt_p.2"
echo "##Samtools index $opt_p.2" >&2
$path_samtools index "$opt_p.2.st.bam"
if [ $? -ne 0 ] || [ ! -s "$opt_p.2.st.bam.bai" ]; then
	echo "Error: Samtools index error: $opt_p.2.st.bam" >&2
	exit 1
fi



### Samtools merge
echo "###Samtools merge $opt_p.1 $opt_p.2"
echo "###Samtools merge $opt_p.1 $opt_p.2" >&2
$path_samtools merge "$opt_p.1.2.st.merge.bam" "$opt_p.1.st.bam" "$opt_p.2.st.bam"
if [ $? -ne 0 ] || [ ! -s "$opt_p.1.2.st.merge.bam" ]; then
	echo "Error: Samtools merge error: $opt_p.1.2.st.merge.bam" >&2
	exit 1
else
	if [ $tempdel -ne 0 ]; then
		rm "$opt_p.1.st.bam" "$opt_p.1.st.bam.bai" "$opt_p.2.st.bam" "$opt_p.2.st.bam.bai"
	fi
fi



### Samtools index
echo "###Samtools index $opt_p.1.2.st.merge.bam"
echo "###Samtools index $opt_p.1.2.st.merge.bam" >&2
$path_samtools index "$opt_p.1.2.st.merge.bam"
if [ $? -ne 0 ] || [ ! -s "$opt_p.1.2.st.merge.bam.bai" ]; then
	echo "Error: Samtools index error: $opt_p.1.2.st.merge.bam" >&2
	exit 1
fi



### picardrmdup
echo "###picardrmdup $opt_p.1.2.st.merge.bam"
echo "###picardrmdup $opt_p.1.2.st.merge.bam" >&2
$path_picardrmdup "$opt_p.1.2.st.merge.bam"
if [ $? -ne 0 ] || [ ! -s "$opt_p.1.2.st.merge.rmdup.bam" ]; then
	echo "Error: Samtools index error: $opt_p.1.2.st.merge.rmdup.bam" >&2;
	exit 1;
else
	if [ $tempdel -ne 0 ]; then
		rm "$opt_p.1.2.st.merge.bam" "$opt_p.1.2.st.merge.bam.bai";
	fi
fi



exit 0;
