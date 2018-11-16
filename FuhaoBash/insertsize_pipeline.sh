#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20160128

Requirements:
	linux
	pcrdup
	bamutils.pcrdup.insertsize.cluster.pl
	perl

Descriptions:
	Estimate inset size, redundancy, and position on chromosome

Options:
  -h    Print this help message
  -i    Input BAM
  -p    output.prefix
  -msd  Maximum position SD from pcrdup
  -isd  Maximum insert SD from pcrdup
  -sb   Insert sizebin for grouping
  -rb   Redundancy sizebin for grouping
  -tb   Position sizebin

Example:
  $0 -i ./my.bam -p output

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
opt_i=''
opt_p=''
maxsd=0
maxinsertsd=0
insertsizebin=0
redundancybin=0
positionbin=0
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -p) opt_p=$2;shift 2;;
    -msd) maxsd=$2;shift 2;;
    -isd) maxinsertsd=$2;shift 2;;
    -sb) insertsizebin=$2;shift 2;;
    -rb) redundancybin=$2;shift 2;;
    -tb) positionbin=$2;shift 2;;
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
str='ni,hai,a'
a=(`SplitString ',' $str`)
echo ${#a[@]} ${a[0]} ${a[1]} ${a[2]}
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
if [ $(CmdExists 'pcrdup') -ne 0 ]; then
	echo "Error: script 'pcrdup' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'bamutils.pcrdup.insertsize.cluster.pl') -ne 0 ]; then
	echo "Error: script 'bamutils.pcrdup.insertsize.cluster.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'sizebin') -ne 0 ]; then
	echo "Error: script 'sizebin' is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################
if [[ "$maxsd" =~ ^[0-9]$ ]]; then
	if [ $maxsd -eq 0 ]; then
		maxsd=20
	elif [ $maxsd -gt 0 ]; then
		echo "*** -msd re-defined to $maxsd"
	else
		echo "***Error: invalid -msd: $maxsd" >&2
	fi
else
	echo "***Error: invalid -msd: $maxsd" >&2
fi
if [[ "$maxinsertsd" =~ ^[0-9]$ ]]; then
	if [ $maxinsertsd -eq 0 ]; then
		maxinsertsd=20
	elif [ $maxinsertsd -gt 0 ]; then
		echo "*** -isd re-defined to $maxinsertsd"
	else
		echo "***Error: invalid -isd: $maxinsertsd" >&2
	fi
else
	echo "***Error: invalid -isd: $maxinsertsd" >&2
fi

if [[ "$insertsizebin" =~ ^[0-9]$ ]]; then
	if [ $insertsizebin -eq 0 ]; then
		insertsizebin=2000
	elif [ $insertsizebin -gt 0 ]; then
		echo "*** -sb re-defined to $insertsizebin"
	else
		echo "***Error: invalid -sb: $insertsizebin" >&2
	fi
else
	echo "***Error: invalid -sb: $insertsizebin" >&2
fi

if [[ "$redundancybin" =~ ^[0-9]$ ]]; then
	if [ $redundancybin -eq 0 ]; then
		redundancybin=5
	elif [ $redundancybin -gt 0 ]; then
		echo "*** -rb re-defined to $redundancybin"
	else
		echo "***Error: invalid -rb: $redundancybin" >&2
	fi
else
	echo "***Error: invalid -rb: $redundancybin" >&2
fi

if [[ "$positionbin" =~ ^[0-9]$ ]]; then
	if [ $positionbin -eq 0 ]; then
		positionbin=10000000
	elif [ $positionbin -gt 0 ]; then
		echo "*** -tb re-defined to $positionbin"
	else
		echo "***Error: invalid -tb: $positionbin" >&2
	fi
else
	echo "***Error: invalid -tb: $positionbin" >&2
fi


#################### Input and Output ###############################
if [ -z "$opt_i" ] || [ ! -s $opt_i ]; then
	echo "Error: invalid input bam" >&2
	exit 1
fi
if [ -z "$opt_p" ]; then
	bamname=${opt_i##*/}
	opt_p==${bamname%.*}
fi



#################### Main ###########################################
echo "##### Input file: $opt_i"
echo "##### Output prefix: $opt_p"
echo "#####	maxsd -msd: $maxsd"
echo "#####	maxinsertsd -isd: $maxinsertsd"
echo "#####	insert sizebin -sb: $insertsizebin"
echo "#####	redundancy sizebin -rb: $redundancybin"
echo "#####	position sizebin -tb: $positionbin"



#Step1: get insert size
echo -e "\n\n\nStep1: get insert size using pcrdup"
pcrdup -counts $opt_p.insertsize.count $opt_i
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count ]; then
	echo "***Error: pcrdup running or output error" >&2
	exit 1
fi



#Step2: deduplicate insert size
echo -e "\n\n\nStep2: deduplicate insert size"
bamutils.pcrdup.insertsize.cluster.pl $opt_p.insertsize.count $maxsd $maxinsertsd > $opt_p.insertsize.count.rmdup
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup ]; then
	echo "***Error: bamutils.pcrdup.insertsize.cluster.pl running or output error" >&2
	exit 1
fi



#Step3: Get insert size
echo -e "\n\n\nStep3: get insert size"
echo "***Total unique insert: "
tail -n +6 $opt_p.insertsize.count.rmdup | wc -l
echo -e "\n***Insert array and counts"
tail -n +6 $opt_p.insertsize.count.rmdup | perl -ne 'chomp;@arr=split(/\t/);print "$arr[2]\t$arr[3]\n";' > $opt_p.insertsize.count.rmdup.insertcount
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.insertcount ]; then
	echo "***Error: substract insert count error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.insertcount"
fi
echo -e "\nEstimate all per pos"
sizebin $opt_p.insertsize.count.rmdup.insertcount $insertsizebin > $opt_p.insertsize.count.rmdup.insertcount.all.sizebin$insertsizebin
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.insertcount.all.sizebin$insertsizebin ]; then
	echo "***Error: sizebin ALL running error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.insertcount.all.sizebin$insertsizebin"
fi
echo -e "\nEstimate unique per pos"
tail -n +6  $opt_p.insertsize.count.rmdup | cut -f 3 > $opt_p.insertsize.count.rmdup.unique
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.unique ]; then
	echo "***Error: substract UNIQUE running error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.unique"
fi
sizebin $opt_p.insertsize.count.rmdup.unique $insertsizebin > $opt_p.insertsize.count.rmdup.unique.sizebin$insertsizebin
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.unique.sizebin$insertsizebin ]; then
	echo "***Error: sizebin UNIQUE error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.unique.sizebin$insertsizebin"
fi



#Step4: redundancy
echo -e "\n\n\nStep4: Estimate redundancy"
tail -n +6 $opt_p.insertsize.count.rmdup | cut -f 4 > $opt_p.insertsize.count.rmdup.redundancy
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.redundancy ]; then
	echo "***Error: redundancy count error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.redundancy"
fi

cat $opt_p.insertsize.count.rmdup.redundancy | perl -ne 'chomp;$sum+=$_;$count++;END {print "***AverageRedundancy: ". $sum/$count ."\n";}'

echo "***Redundancy Bin: $redundancybin"
sizebin $opt_p.insertsize.count.rmdup.redundancy $redundancybin > $opt_p.insertsize.count.rmdup.redundancy.bin$redundancybin
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.redundancy.bin$redundancybin ]; then
	echo "***Error: redundancy sizebin error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.redundancy.bin$redundancybin"
fi



##Step6: Chromosome distribution
echo -e "\n\n\nStep5: get position distribution"
tail -n +6 $opt_p.insertsize.count.rmdup | cut -f 2 > $opt_p.insertsize.count.rmdup.posn
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.posn ]; then
	echo "***Error: substract mapping position error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.posn"
fi
echo "***Position Bin $positionbin"
sizebin $opt_p.insertsize.count.rmdup.posn $positionbin > $opt_p.insertsize.count.rmdup.posn.bin$positionbin
if [ $? -ne 0 ] || [ ! -s $opt_p.insertsize.count.rmdup.posn.bin$positionbin ]; then
	echo "***Error: substract mapping position error" >&2
	exit 1
else
	echo "***Data stored in file: $opt_p.insertsize.count.rmdup.posn.bin$positionbin"
fi

echo "ALL DONE"

exit 0
