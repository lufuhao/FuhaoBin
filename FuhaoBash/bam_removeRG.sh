#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi
machtype=${machtype%% *}


################# help message ######################################
help() {
cat<<HELP

$0 --- Remove @RG

Version: 20150518

Requirements:
	samtools
	perl

Descriptions:
	Remove @RG in header and RG:Z:xx in alignments from SAM/BAM
	*Auto-detect SAM/BAM based on extension .sam or .bam
	*Input in SAM/BAM, output in BAM
	*Output default: inbam_basename.noRG.bam

Options:
  -h    Print this help message
  -i    CONFIG file
  -t    Number of threads, default: 1
  -s    Not run simulation
  -a    Not run assembly

Example:
  $0 -i in.bam/in.sam -o out.bam

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
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH



#################### Initializing ###################################
inbam=''
outbam=''
debug=0

#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) inbam=$2;shift 2;;
    -o) outbam=$2;shift 2;;
    -d) debug=1; shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done

#################### Subfuctions ####################################
###Detect command existence
function CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0
  else
    echo 1
  fi
}



#################### Command test ###################################
if [ $(CmdExists 'samtools') -eq 1 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'samtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################
if [ -z "$inbam" ] || [ ! -s $inbam ]; then
	echo "Error: invalid BAM input: $inbam" >&2
	exit 1
fi
inbm_name=${inbam##*/}
inbm_base=${inbm_name%.*}
inbm_exts=${inbm_name##*.}
if [[ -z "$outbam" ]]; then
	outbam="$inbm_base.noRG.bam"
	if [ $debug -eq 1 ]; then
		echo "Info: Output BAM to $outbam"
	fi
fi
if [ -e $outbam ]; then
	echo "Error: invalid or existing BAM output: $outbam" >&2
	exit 1
fi



#################### Main ###########################################
if [[ $inbm_exts =~ [sS][aA][mM] ]]; then
	if [ $debug -eq 1 ]; then
		echo "Info: BAM input in SAM format"
	fi
	samtools view -S -h $inbam | perl -ne 'chomp; $line=$_; if ($line=~/^\@/) { print $line."\n" unless ($line=~/^\@RG/);}else{if ($line=~/\s+RG:Z\S+/) {$line=~s/\s+RG:Z\S+//g;} print $line."\n";}' | samtools view -b -S -h - > $outbam
elif [[ $inbm_exts =~ [bB][aA][mM] ]]; then
	if [ $debug -eq 1 ]; then
		echo "Info: BAM input in BAM format"
	fi
	samtools view -h $inbam | perl -ne 'chomp; $line=$_; if ($line=~/^\@/) { print $line."\n" unless ($line=~/^\@RG/);}else{if ($line=~/\s+RG:Z\S+/) {$line=~s/\s+RG:Z\S+//g;} print $line."\n";}' | samtools view -b -S -h - > $outbam
	if [ $? -ne 0 ]; then
		echo "Error: samtools parsing failed" >&2
		exit 1
	fi
else
	echo "Error: unknown BAMinput extension: $inbm_exts" >&2
	exit 1
fi

if [ -e $outbam.bai ]; then
	rm $outbam.bai > /dev/null 2>&1
fi
samtools index $outbam
if [ $? -ne 0 ] || [ ! -s $outbam.bai ]; then
	echo "Error: samtools index error: $outbam" >&2
	exit 1
fi

exit 0

