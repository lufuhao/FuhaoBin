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

################# help message ######################################
help() {
cat<<HELP

echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"

$0 --- running exonerate in server mode

Version: 20170707

Requirements:
  Linux: perl, grep
  Script: 
      fasta_splitter_by_numseq.pl
      exonerate_gff2_to_gff3.pl
  Exonerate

Descriptions:
  Split fasta: NUM [-n] seq each
  Map sequences to reference in server mode
  Collect Exonerate out [GTF/gff2 format]
  Transform gtf to gff3

Options:
  -h    Print this help message
  -i    Query ESTs in fasta
  -g    Genome seqences in fasta
  -o    Final output (Default: \$prefix.EXONERATE_OUT)
  -n    Seqeunce number to split for each fasta file (Default: 5000)
  -e    Exenonerate --percent option (Default: 70)
  -p    Output prefix (Default: "MyExonerate")
  -m    Memory limit (MB) (Default: 1024)
  -b    -bestn [INT] for exonerate [10]
  -r    Server port number (Default: 12886)
  -s    Skip server setup
          * Will diable -g -m options
          fasta2esd -s FALSE -f \$opt_g -o \$opt_p.esd
          esd2esi \$opt_p.esd \$opt_p.esi --memorylimit \$opt_m
          nohup exonerate-server \$opt_p.esi --port \$opt_r  &
  -x    close server
  -k    Keep temporary file

Example:
  $0 -i cds.fasta -g genome.fa -o exonerateout -n 500 -m 20000 -r 12886 -s

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
echo -e "\n#########\n${ProgramName} Program initializing ...\n###########\n"


#################### Initializing ###################################
opt_i=""
opt_g=""
opt_n=5000
opt_p="MyExonerate"
opt_m=1024
opt_r=12886
opt_o=0
opt_s=0
opt_e=70
opt_x=0
bestn=10
opt_k=0;
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -g) opt_g=$2;shift 2;;
    -o) opt_o=$2;shift 2;;
    -n) opt_n=$2;shift 2;;
    -e) opt_e=$2;shift 2;;
    -p) opt_p=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    -r) opt_r=$2;shift 2;;
    -b) bestn=$2;shift 2;;
    -k) opt_k=1;shift;; 
    -s) opt_s=1;shift;;
    -x) opt_x=1;shift;;
    --) shift;break;;
    -*) echo "${ProgramName}Error: no such option $1. -h for help" > /dev/stderr;exit 1;;
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
if [ $(CmdExists 'perl') -ne 0 ]; then
	echo "${ProgramName} Error: CMD 'perl' in PROGRAM 'PERL' not found" >&2 
	exit 127
fi
if [ $(CmdExists 'grep') -ne 0 ]; then
	echo "${ProgramName} Error: CMD 'grep' not found." >&2 
	exit 127
fi
if [ $(CmdExists 'fasta_splitter_by_numseq.pl') -ne 0 ]; then
	echo "${ProgramName} Error: script 'fasta_splitter_by_numseq.pl' not found." >&2 
	exit 127
fi
if [ $(CmdExists 'exonerate') -ne 0 ] ; then
	echo "${ProgramName} Error: CMD 'exonerate' in PROGRAM 'Exonerate' not found." >&2 
	exit 127
fi
if [ $opt_s -eq 0 ]; then
	if [ $(CmdExists 'fasta2esd') -ne 0 ] ; then
		echo "${ProgramName} Error: CMD 'fasta2esd' in PROGRAM 'Exonerate' not found." >&2 
		exit 127
	fi
	if [ $(CmdExists 'esd2esi') -ne 0 ] ; then
		echo "${ProgramName} Error: CMD 'esd2esi' in PROGRAM 'Exonerate' not found." >&2 
		exit 127
	fi
	if [ $(CmdExists 'exonerate-server') -ne 0 ] ; then
		echo "${ProgramName} Error: CMD 'exonerate-server' in PROGRAM 'exonerate-server' not found." >&2 
		exit 127
	fi
fi
if [ $(CmdExists 'exonerate_gff2_to_gff3.pl') -ne 0 ] ; then
	echo "${ProgramName} Error: SCRIPT 'exonerate_gff2_to_gff3.pl' not found." >&2 
	exit 127
fi




#################### Defaults #######################################
if [ -z "$opt_i" ]; then
	echo "${ProgramName}Error: input ESTs sequence file not specified" >&2;
	exit 100;
fi
opt_i=$(echo $(cd $(dirname "$opt_i"); pwd)/$(basename "$opt_i"))
if [ ! -s "$opt_i" ]; then
	echo "${ProgramName}Error: input ESTs sequence file not found" >&2;
	exit 100;
fi
if [ $opt_s -eq 0 ]; then
	if [ -z "$opt_g" ] || [ ! -s "$opt_g" ]; then
		echo "${ProgramName}Error: input genome sequence file not found" >&2;
		exit 100;
	fi
	if [[ "$opt_m" =~ ^[0-9]+$ ]] && (( $opt_m > 0 )); then
		echo "${ProgramName}Info: memory limit accepted: $opt_m"
	else
		echo "${ProgramName}Error: invalid memory limit specified by '-m'" >&2;
		exit 100;
	fi
fi
if [[ "$opt_n" =~ ^[0-9]+$ ]] && (( $opt_n > 0 )); then
	echo "${ProgramName}Info: split number accepted: $opt_n"
else
	echo "${ProgramName}Error: invalid split number specified by '-n'" >&2;
	exit 100;
fi

if [[ "$opt_r" =~ ^[0-9]+$ ]] && (( $opt_r > 0 )); then
	echo "${ProgramName}Info: port number accepted: $opt_r"
else
	echo "${ProgramName}Error: invalid port number specified by '-r'" >&2;
	exit 100;
fi
if [ -z "$opt_p" ]; then
	echo "${ProgramName}Error: invalid output prefix specified by '-p'" >&2;
	exit 100;
fi

if [ "$opt_o" = "0" ]; then
	opt_o="${opt_p}.exonerateout"
elif [ -z "$opt_o" ]; then
	opt_o="${opt_p}.exonerateout"
	echo "${ProgramName}Error: invalid final output specified by '-o'" >&2;
	echo "                     using default -o : $opt_o" >&2;
fi
if [ -e "$opt_o" ]; then
	echo "${ProgramName}Error: final output existing: $opt_o" >&2;
	exit 100;
fi


#################### Input and Output ###############################




#################### Main ###########################################
if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "${ProgramName}Error: sort error" >&2
	exit 1
fi

rundir=$PWD

echo "### PWD: $rundir"
if [ $opt_s -eq 0 ]; then
	echo "### 1. Establish Exonerate server: port $opt_r"
	fasta2esd -s FALSE -f $opt_g -o $opt_p.esd
	if [ $? -ne 0 ] || [ ! -s "$opt_p.esd" ]; then
		echo "${ProgramName}Error: building database error" >&2
		exit 1
	fi
	esd2esi $opt_p.esd $opt_p.esi --memorylimit $opt_m
	if [ $? -ne 0 ] || [ ! -s "$opt_p.esi" ]; then
		echo "${ProgramName}Error: indexing database error" >&2
		exit 1
	fi
	#exonerate-server $opt_p.esi --port $opt_r &
	nohup exonerate-server $opt_p.esi --port $opt_r  &
	if [ $? -ne 0 ]; then
		echo "${ProgramName}Error: build exonerate server error" >&2
		exit 1
	fi
fi

echo "### 2. Spilt fasta: $opt_n seq each file"
if [ -d "$rundir/split" ]; then
	echo "${ProgramName} Warnings: temporary folder exists: $rundir/split" >&2
	echo "${ProgramName} Warnings: using current fasta splits..." >&2
else
	mkdir -p $rundir/split
	cd $rundir/split
	fasta_splitter_by_numseq.pl $opt_i $opt_n $opt_p
fi



cd $rundir/
if [ -d "$rundir/exonerateout" ]; then
	echo "Error: temporary folder exists: $rundir/exonerateout" >&2
	exit 1
else
	mkdir -p $rundir/exonerateout
fi
cd $rundir/exonerateout
for myfasta in `ls $rundir/split/${opt_p}*.fa`; do
	echo "### Mapping fasta: $myfasta"
	fabasename=${myfasta##*/}
	exonerate --model est2genome $myfasta localhost:$opt_r --percent $opt_e --score 100 --showvulgar yes --bestn $bestn --minintron 20 --softmaskquery no --softmasktarget no --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --showalignment no  --showtargetgff yes --geneseed 250  > $rundir/exonerateout/$fabasename.EXONERATE_OUT
	if [ $? -ne 0 ] || [ ! -s "$rundir/exonerateout/$fabasename.EXONERATE_OUT" ]; then
		echo "${ProgramName}Warnings: Exonerate might fail on $myfasta" >&2
	fi
	
done



if [ $opt_x -eq 1 ]; then
	pidnum=$(ps -ef | grep 'exonerate-server' | grep "port $opt_r" | perl -ne 'BEGIN{@pros=();} chomp; @arr=split(/\s+/); if ($arr[1] =~ /^\d+$/) {push (@pros, $arr[1]);} END {if (scalar(@pros)==1) {print $pros[0], "\n ";}}')
	if [ ${#pidnum[@]} -eq 1 ] && [[ "${pidnum[0]}" =~ [0-9]+ ]]; then
		kill ${pidnum[0]}
	fi
fi


cd $rundir/
cat $rundir/exonerateout/*EXONERATE_OUT > $opt_o
if [ $? -ne 0 ] || [ ! -s "$opt_o" ]; then
	echo "${ProgramName}Error: Output error: $opt_o" >&2
	exit 100;
else
	exonerate_gff2_to_gff3.pl $opt_o > $opt_o.gff3
	if [ $? -ne 0 ] || [ ! -s "$opt_o" ]; then
		echo "${ProgramName}Error: exonerateout to GFF3 conversion error: $opt_o.gff3" >&2
		exit 100
	fi
fi

if [ $opt_k -ne 0 ]; then
	echo "${ProgramName}info: Cleaning $rundir/exonerateout $rundir/split"
	rm -rf $rundir/exonerateout $rundir/split > dev/null 2>&1
fi

echo "### JOB done ###"
exit 0;
