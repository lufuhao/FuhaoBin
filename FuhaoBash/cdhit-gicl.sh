#!/bin/bash
RootDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)


################# help message ######################################
help () {
cat<<HELP

cdhit-gicl -i input.fa [OPTIONS]

Version: 20181115

Descriptions:
  run cdhit-GICL-cdhit pipeline
  
Require:
	cd-hit-est
	gicl.pl
	ace2fasta.pl
	cdbyank

Options:
  -h    Print this help message
CD-HIT
  -i    Redundant fasta file
  -o    output fasta file
  -c    Sequence identity threshold, default 0.9
  -n    Word_length, default 10
  -t    Number of threads, default 1; with 0, all CPUs will be used
  -r    1 or 0, default 1, by default do both +/+ & +/- alignments
        if set to 0, only +/+ strand alignment
  -m    Memory limit (in MB) for the program, default 800; 0 for unlimitted
GICL
  -u    Use the specified number of CPUs on local machine; default 1
  -l    Miminum overlap length (default 40)
  -p    Minimum percent identity for overlaps (default 94)
  -v    Maximum length of unmatched overhangs (default 30)
  -k    Lower-case masking in <fasta_db> sequences
  -s    Number of sequences in a clustering search slice (default 1000)

Example:
  cdhit-gicl -i input.fa -o output.fa -c 0.95 -n 8 -t 0 -r 1 -m 20000 -u 10 -l 50 -p 99 -v 100 -s 2000
20150320:  cdhit-gicl -i input.fa -o output.fa -c 0.95 -n 8 -t 0 -r 1 -m 20000 -u 10 -l 200 -p 95 -v 30 -s 2000

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
#################### Defaults #######################################
echo -e "\n######################\ncdhit-GICL initializing ...\n######################\n"



###Defaults##########################################################
opt_c=0.9
opt_n=10
opt_t=0
opt_r=1
opt_m=800
opt_u=1
opt_l=40
opt_p=94
opt_v=30
opt_k=1
opt_s=1000
opt_o=''

#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -o) opt_o=$2;shift 2;;
    -c) opt_c=$2;shift 2;;
    -n) opt_n=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -r) opt_r=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    -u) opt_u=$2;shift 2;;
    -l) opt_l=$2;shift 2;;
    -p) opt_p=$2;shift 2;;
    -v) opt_v=$2;shift 2;;
    -s) opt_s=$2;shift 2;;
    -k) opt_k=0;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
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



#################### Command test ###################################
if [ $(CmdExists 'cd-hit-est') -eq 1 ]; then
	echo "Error: CMD cd-hit-est in PROGRAM 'CD-HIT'  is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'gicl.pl') -eq 1 ]; then
	echo "Error: script 'gicl.pl' in PROGRAM 'GICL'  is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'ace2fasta.pl') -eq 1 ]; then
	echo "Error: script ''ace2fasta.pl'' in PROGRAM 'GICL'  is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'cdbyank') -eq 1 ]; then
	echo "Error: CMD 'cdbyank' in PROGRAM 'CDBtools'  is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################
if [ -e $opt_i ]; then
  if [ ! -s $opt_i ]; then
    echo "Error: empty input" >&2
    exit 1
  fi
else
  echo "Error: invalid input, file not exists" >&2
  exit 1
fi
filename=${opt_i##*/}
fileext=${filename##*.}
filebasename=${filename%.*}
filetissue=${filename%%_*}

if [ -z "$opt_o" ];then
  opt_o="${filetissue}_RemoveRedundancy.fa"
fi



#################### Main ###########################################
RunDir=$(pwd)
echo "Output path: $RunDir"
###1.cdhit
if [ ! -d $RunDir/1.$filetissue.cdhit ]; then
  mkdir -p $RunDir/1.$filetissue.cdhit
fi
cd $RunDir/1.$filetissue.cdhit
echo "CMD: cd-hit-est -i $opt_i -c $opt_c -n $opt_n -T $opt_t -r $opt_r -M $opt_m -o $filetissue.cdhit.fa"
cd-hit-est -i $opt_i -c $opt_c -n $opt_n -T $opt_t -r $opt_r -M $opt_m -gap -2 -o $filetissue.cdhit.fa

if [ -s $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa ]; then
  if [ ! -d $RunDir/2.$filetissue.gicl ]; then
    mkdir -p $RunDir/2.$filetissue.gicl
  fi
else
  echo "Error: CDHIT output error"
  exit 1
fi


###2.GICL
cd $RunDir/2.$filetissue.gicl
if [ $opt_k -eq 1 ];then
  echo "CMD: gicl.pl -c $opt_u -in $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa -l $opt_l -p $opt_p -v $opt_v -mk -n $opt_s"
  gicl.pl -c $opt_u -in $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa -l $opt_l -p $opt_p -v $opt_v -mk -n $opt_s
elif [ $opt_k -eq 0 ];then
  echo "CMD: gicl.pl -c $opt_u -in $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa -l $opt_l -p $opt_p -v $opt_v -n $opt_s"
  gicl.pl -c $opt_u -in $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa -l $opt_l -p $opt_p -v $opt_v -n $opt_s
fi
if [ ! -s $RunDir/2.$filetissue.gicl/*.ace ]; then
  echo "ace2fasta.pl output error ..."
  exit 1
fi
ace2fasta.pl -o contigs.fa *.ace
if [ ! -s $RunDir/2.$filetissue.gicl/contigs.fa ]; then
  echo "File contigs.fa does not exist"
  exit 1
fi

if [ ! -s $RunDir/2.$filetissue.gicl/*.singletons ]; then
  echo "File contigs.fa does not exist"
  exit 1
fi
cat *.singletons | cdbyank $RunDir/1.$filetissue.cdhit/$filetissue.cdhit.fa.cidx > singletons.fa
if [ ! -s $RunDir/2.$filetissue.gicl/singletons.fa ]; then
  echo "cdbyank output singletons.fa error"
  exit 1
fi
cat contigs.fa singletons.fa > $RunDir/2.$filetissue.gicl/$filetissue.gicl.fa
if [ -s $RunDir/2.$filetissue.gicl/$filetissue.gicl.fa ]; then
  if [ ! -d $RunDir/3.$filetissue.cdhit ]; then
    mkdir -p $RunDir/3.$filetissue.cdhit
  fi
else
  echo "Error: GICL output error"
  exit 1
fi


###3.cdhit
cd $RunDir/3.$filetissue.cdhit
echo "CMD: cd-hit-est -i $RunDir/2.$filetissue.gicl/$filetissue.gicl.fa -c $opt_c -n $opt_n -T $opt_t -r $opt_r -M $opt_m -gap -2 -o 3.$filetissue.cdhit.fa"
cd-hit-est -i $RunDir/2.$filetissue.gicl/$filetissue.gicl.fa -c $opt_c -n $opt_n -T $opt_t -r $opt_r -M $opt_m -gap -2 -o 3.$filetissue.cdhit.fa
if [ ! -s $RunDir/3.$filetissue.cdhit/3.$filetissue.cdhit.fa ]; then
  echo "File 3.$filetissue.cdhit.fa does not exist"
  exit 1
fi
mv 3.$filetissue.cdhit.fa $RunDir/$opt_o


###4.output
cd $RunDir/
#rm -rf $RunDir/1.$filetissue.cdhit $RunDir/2.$filetissue.gicl $RunDir/3.$filetissue.cdhit
echo "CDHIT-gicl finished"


exit 0
