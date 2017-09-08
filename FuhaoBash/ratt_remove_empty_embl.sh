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

Version: 20170324

Requirements:

Descriptions:
	Clean empty EMBL file after RATT

Example:
  $0 EMBL_dir empty_dir

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
#echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

embldir=$1
newdir=$2

if  [ -z "$newdir" ]; then
	newdir="empty"
fi

if [ ! -d $newdir ]; then
	mkdir -p $newdir
fi


for emblfile in `find $embldir/ -name *.final.embl -type f`; do
#	echo $emblfile
	ftlines=0
	ftlines=$(grep ^'FT' $emblfile | wc -l)
#	if [ $ftlines -gt 0 ]; then
#		echo "$ftlines lines"
#	fi
	if [ $ftlines -lt 1 ]; then
		mv ${emblfile%.final.embl}* $newdir
	fi
done

