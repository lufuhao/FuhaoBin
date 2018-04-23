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
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"


################# help message ######################################
help() {
cat<<HELP

$0 --- Sort GFF

Version: 20170612

Requirements:
	Linux: grep, sort
	SAMtools: tabix, bgzip

Descriptions:
	Sort, BGzip, and index GFF file

Options:
  -h    Print this help message
  -i    GFF3 input
  -o    [Opt] GFF.gz output; Default: inputbase.sorted.gff.gz

Example:
  $0 -i my.gff3 -o my.gff3.gz

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

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) gffin=$2;shift 2;;
    -o) gffout=$2;shift 2;;
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

    echo 1
  fi
}



#################### Command test ###################################
if [ $(CmdExists 'grep') -eq 1 ]; then
	echo "Error: CMD 'grep' on Linux is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'sort') -eq 1 ]; then
	echo "Error: CMD 'sort' on Linux is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'bgzip') -eq 1 ]; then
	echo "Error: CMD 'bgzip' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127;
fi
if [ $(CmdExists 'tabix') -eq 1 ]; then
	echo "Error: CMD 'tabix' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127;
fi



#################### Defaults #######################################




#################### Input and Output ###############################
if [ -z "$gffin" ]; then
	echo "GFFSORT_Error: need to specify GFF input" >&2
	exit 100
elif [ ! -s "$gffin" ]; then
	echo "GFFSORT_Error: invalid GFF input" >&2
	exit 100
fi
if [ -z "$gffout" ]; then
#	echo "GFFSORT_Info: GFF output not specified, use default instead"
	gffin_name=${gffin##*/}
	gffin_base=${gffin%.*}
	gffout="$gffin_base.sorted.gff.gz"
	echo "GFFSORT_Info: using default GFF output: $gffout"
fi
if [ -e $gffout ]; then
	echo "GFFSORT_Warnings: existing GFF output" >&2
	echo "GFFSORT_Warnings: Rename existing $gffout to $gffout.bak" >&2
	mv -i $gffout $gffout.bak
fi



#################### Main ###########################################
echo -e "\n### GFFSORT_Info: Sorting"
(grep ^"#" $gffin; grep -v ^"#" $gffin | sort -k1,1 -k4,4n) | bgzip > $gffout;
if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	exit 100
fi
tabix -p gff $gffout;
echo -e "\n### GFFSORT_Info: tabix index GFF"
if [ $? -ne 0 ] || [ ! -s $gffout.tbi ]; then
	echo "GFFSORT_Error: invalid GFF input" >&2
	exit 100
fi
#tabix sorted.gff.gz chr1:10,000,000-20,000,000; 
echo -e "\n### GFFSORT_Info: Finished\n";



exit 0;
