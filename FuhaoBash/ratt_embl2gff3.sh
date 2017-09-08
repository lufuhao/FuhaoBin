#!/bin/bash
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

$0 --- convert RATT output EMBL to GFF3

Version: 20170620

Requirements:
	embl2gff3.pl

Descriptions:
	Input: RATT EMBL dir [default: ./transfered ]
	Output: GFF3 dir [default: ./transferedgff3 ]

Options:
  -h    Print this help message
  -i    embl_dir
  -o    out_dir

Example:
  $0 -i embl_dir -o out_dir

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
embldir=''
gffdir=''
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) embldir=$2;shift 2;;
    -o) gffdir=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" >&2;exit 1;;
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
if [ $(CmdExists 'embl2gff3.pl') -ne 0 ]; then
	echo "Error: script 'embl2gff3.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi




#################### Inout and Output ###############################
if [ -z "$embldir" ]; then
	embldir='./transfered'
	echo "Info: using default input dir: $embldir"
elif [ -d "$embldir" ]; then
	echo "Info: using specified input dir: $embldir"
else
	echo "Error: EMBL dir not found\n" >&2
	exit 10;
fi

if [ ! -z "$gffdir" ]; then
	echo "Info: using specified output dir: $gffdir"
else
	gffdir='./transferedgff3'
	echo "Info: using default output dir: $gffdir"
fi

if [ -d "$gffdir" ]; then
	echo "    * Existing output dir"
	echo "    * Cleaning output dir"
	rm -rf "$gffdir"/* > /dev/null 2>&1
else
	echo "    * Create output dir"
	mkdir -p "$gffdir"
fi




for emblfile in `find $embldir/ -name *.final.embl`; do
	gff3file=${emblfile##*/}
	echo -e "\n\n\n$gff3file"
	embl2gff3.pl $emblfile $gffdir/$gff3file.gff3
	if [ $? -ne 0 ]; then
		echo "$emblfile failed" >&2
	fi
done
