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

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20150603

Requirements:
	samtools
	HISAT2
	stringtie
	ballgown

Descriptions:
	xxx

Options:
  -h    Print this help message
  -i    CONFIG file
  -t    Number of threads, default: 1
  -s    Not run simulation
  -a    Not run assembly

Example:
  $0 -i ./chr1.fa -t 10

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
opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) opt_s=1;shift 1;;
    -a) opt_a=1;shift 1;;
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
if [ $(CmdExists 'santools') -eq 0 ]; then
	exit 0
else
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'santools') -eq 1 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################
outputpath=
outprefix="ParD27G1"
readgroup="ParD27G1"
read1=
read2=
threads=
hisat2index=
genome=
#################### Main ###########################################
#./hisat2-2.0.0-beta/hisat2-build -f ucsc.hg19.fasta --ss splicesites.txt --exon exons.txt -p 7 ./ucsc.hg19
## 添加--ss和--exon选项后，需要很大的内存，build 人基因组的话需要200G RAM，如果没有这么大内存，不要添加这两个选项，但要在后续运行hisat时添加 --known-splicesite-infile选项
#hisat2-build -f ucsc.hg19.fasta -p 7 ./uscs.hg19 ##大概需要一小时二十分钟
#build index的过程和bowtie2很像。
#-q指定fastq格式
#hisat2 -q -x ./ucsc.hg19 -1 reads_1.fastq -2 reads_2.fastq -S alns.sam -t
## hisat2 -q -x ./ucsc.hg19 -1 reads_1.fastq -2 reads_2.fastq -S alns.sam --known-splicesite-infile splicesites.txt -t
echo "Step1: HISAT2"
echo "Readgroup: $readgroup"
echo "OutBAM: $outprefix" >&2
#without splicing sites and exons
hisat2-build -p $threads $genome $hisat2index
#HISAT2
hisat2 -q -p $threads --fr \
--novel-splicesite-outfile $outprefix.splicesites2.new \
--rg-id "$readgroup" --rg "SM:$readgroup" --rg "PL:Illumina" \
-x Aet3_hisat2 -1 $read1 -2 $read2 | samtools view -b -h -S -F 4 - > $outputpath/$outprefix.$hisat2index.bam
if [ $? -ne 0 ] || [ ! -s "$outputpath/$outprefix.$hisat2index.bam" ]; then
	echo "GFFSORT_Error: sort error" >&2
	exit 1
}


if [ ! -d $sortdir ]; then
	mkdir -p $sortdir
fi
cd $sortdir
echo "###BAMSORT"
for bamfile in `ls $bamfilepath/*.bam`; do
	bamname=${bamfile##*/}
	bambase=${bamname%.bam}
	echo $bamname
	echo $bamname >&2
	samtools sort $bamfile $bambase.sort
	if [ $? -ne 0 ] || [ ! -s "$bambase.sort.bam" ]; then
		echo "Error: bam sort" >&2
		continue
	fi
	samtools index $bambase.sort.bam
	if [ $? -ne 0 ] || [ ! -s "$bambase.sort.bam" ]; then
		echo "Error: bam index" >&2
		continue
	fi
	echo "$bambase.sort.bam"
	echo "$bambase.sort.bam" >&2
	samtools flagstat $bambase.sort.bam
done

if [ ! -d "$deduppath" ]; then
	mkdir -p "$deduppath"
fi
cd "$deduppath"
picardrmdup -a $sortdir/*.sort.bam



