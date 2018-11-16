#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi
abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
#NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"

################# help message ######################################
help() {
cat<<HELP

$0 --- extract seq from fasta

Version: 201500901

Requirements:
    CDBtools

Descriptions:
    Extract fasta using cdbfasta and cdbyank

Options:
  -h    Print this help message
  -i    Fasta database
  -s    ID list, comma delimited
  -f    IF list file
  -p    Outpath
  -o    Outname

Example:
  $0 -i ./chr1.fa -s seq1,seq2 -o out.fa

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
fastafile=''
idlist=''
idfile=''
outpath=''
outfa=''
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) fastafile=$2;shift 2;;
    -s) idlist=$2;shift 2;;
    -f) idfile=$2;shift 2;;
    -p) outpath=$2;shift 2;;
    -o) outfa=$2; shift 2;;
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
}
###Usage: array=(`split delimiter string`)
split () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}
#str='ni,hai,a'
#a=(`SplitString ',' $str`)
#echo ${#a[@]} ${a[0]} ${a[1]} ${a[2]}
#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}



#################### Command test ###################################
if [ $(CmdExists 'cdbfasta') -ne 0 ]; then
	echo "Error: CMD 'cdbfasta' in PROGRAM 'CDBtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'cdbyank') -ne 0 ]; then
	echo "Error: CMD 'cdbyank' in PROGRAM 'CDBtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
declare -a idssum=()
declare -a id1=()
declare -a id2=()
if [ -z "$idlist" ]; then
	echo "Info: not using ID list -s"
else
	echo "Info: using ID list -s $idlist"
	id1=($(echo $idlist | perl -ne 'next unless (/^\S+$/); @arr=split(/,/); foreach (@arr) {print $_."\n";}'))
	echo "Info: total ${#id1[@]} IDs: using ID list -s: ${id1[@]}"
fi
if [ -z "$idfile" ] || [ ! -s "$idfile" ]; then
	echo "Info: not using ID file -f"
else
	echo "Info: using ID list -f $idfile"
	id2=($(cat $idfile))
	echo "Info: total ${#id2[@]} IDs: using ID list -s: ${id2[@]}"
fi
idssum=("${id1[@]}" "${id2[@]}")
if [ ${#idssum[@]} -lt 1 ]; then
	echo "Error: empty IDs to extract" >&2
	exit 1
else
	echo "Info: SUM ${#idssum[@]} IDs: ${idssum[@]}"
	cdbyankid=$(echo "${idssum[@]}" | perl -ne 's/^\s+//; s/\s+$//; next unless (/\S+/); s/\s+/\\n/g;print;')
	echo "Test: cdbyankid: $cdbyankid"
fi
#################### Defaults #######################################
output=''



#################### Input and Output ###############################
if [ -z "$fastafile" ] || [ ! -s "$fastafile" ]; then
	echo "Error: invalid fasta input" >&2
	exit 1
fi
fastaname=${fastafile##*/}

if [ ! -s "$fastafile.cdbz" ]; then
	echo "Info: use cdbfasta to index fasta $fastafile"
	cdbfasta $fastafile -z $fastafile.cdbz
	if [ -s "$fastafile.cdbz" ]; then
		echo "Info: successfully use cdbfasta to index fasta $fastafile"
	else
		echo "Error: Failed to use cdbfasta to index fasta $fastafile" >&2
		exit 1
	fi
else
	echo "Info: use cdbfasta to index fasta $fastafile"
fi
if [ ! -s "$fastafile.cdbz.cidx" ]; then
	echo "Error: fasta index not exist : $fastafile.cdbz.cidx" >&2
	exit 1
fi
if [ -z "$outpath" ]; then
	echo "Info: out path (-o) not specified"
	if [ -z "$outfa" ]; then
		echo "Error: both output path and filename not specified"
		exit 1
	else
		output=$outfa
	fi
else
	if [ ! -d "$outpath" ]; then
		echo "Info: create Directory: $outpath"
		mkdir -p "$outpath"
	fi
	if [ -z "$outfa" ]; then
		output="$outpath/$fastaname.extract.fa"
	else
		output="$outpath/$outfa"
	fi
fi



#################### Main ###########################################
echo -e "\n\n\n###Summary ###"
echo -e "\tInputDB: $fastafile"
echo -e "\tIDs: $cdbyankid"
echo -e "\tOutput: $output"

### extract fasta
echo -e "$cdbyankid" | cdbyank $fastafile.cdbz.cidx -o $output -w
if [ $? -ne 0 ] || [ ! -s $output ]; then
	echo "cdbyank: sort error" >&2
	exit 1
fi

exit 0
