#!/bin/sh
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

repeatscout-pipe.sh --- NGS sIMulation PipeLinE

Version: 20180412

Descriptions:
  Run RepeatScout to discover tandem repeats

Requirements
  Program: RepeatScout
               build_lmer_table
               RepeatScout
               filter-stage-1.prl
               filter-stage-2.prl
           RepeatMasker

Options:
  -h      Print this help message
  -i      Fasta file
  -t     [1] The number of processors in parallel 
        (only works for batch files or sequences over 50 kb)
  -n     [10] The number of times a sequence must appear
           for filter-stage-2.prl to be reported. 
  -rm      ["RepeatMasker"] Path/to/Repeatmasker
  -c      Run final RepeatMasker
  -p      [fasta basename] Output prefix 



Example:
  repeatscout-pipe.sh -i my.fa -t 1 -n 10 -p MyOut

Output 
  5.OutputPrefix.secondfilter.lib

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
[ "$1" = "-h" ] && help


#################### Defaults #######################################
echo -e "\n######################\nRun RepeatScout initializing ...\n######################\n"


###Defaults
num_cpus=1
thres=10
cmd_repeatmask='RepeatMasker'
run_final_repeatmasker=0
opt_p="MyOut"
#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -t) num_cpus=$2;shift 2;;
    -n) thres=$2;shift 2;;
    -rm) cmd_repeatmask=$2;shift 2;;
    -c) run_final_repeatmasker=1;shift 1;;
    -p) opt_p=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help";exit 1;;
    *) break;;
  esac
done




###Check if necessary commands exist
function CmdExists {
#  command -v $1 >/dev/null 2>&1 || { echo >&2 "I require $1 but it's not installed.  Aborting."; exit 1; }
  local cmd=$1
  if command -v $cmd >/dev/null 2>&1;then
    echo $cmd "  :  "`command -v $cmd` >&2
    exit 0
  else
    echo "Error: require $cmd but it's not installed.  Exiting..." >&2
    exit 100
  fi
}
if [ `CmdExists build_lmer_table` ] || [ `CmdExists RepeatScout` ] || [ `CmdExists filter-stage-1.prl` ] || [ `CmdExists filter-stage-2.prl` ]; then
  echo "### Please install RepeatScout or put the executables into PATH first" >&2
  exit 100
fi
if [ `CmdExists RepeatMasker` ]; then
  echo "### Please install RepeatMasker or put the executables into PATH first" >&2
  exit 100
fi

### Input and output
if [ ! -s $opt_i ]; then
  echo "Fasta input not found ..." >&2
  exit 100
fi

fasta_name=${opt_i##*/}
fasta_base=${fasta_name%.*}
if [ "$opt_p" = "MyOut" ]; then
	opt_p=$fasta_base
fi
echo "Info: input fasta  : $opt_i"
echo "Info: Output prefix: $opt_p"

### Clean previous output
if [ -s $opt_p.1.12base.freq ] || [ -s $opt_p.2.firstrepeats.lib ] || [ -s $opt_p.3.filtered.lib ] || [ -s $opt_p.4.rm1.out ] || [ -s $opt_p.5.secondfilter.lib ] || [ -s $opt_p.6.rm2.out ]; then
  rm $opt_p.1.12base.freq $opt_p.2.firstrepeats.lib $opt_p.3.filtered.lib $opt_p.4.* $opt_p.5.secondfilter.lib $opt_p.6.*
fi


### 1. build frequency table
### RepeatScout: build_lmer_table 
build_lmer_table -sequence $opt_i -freq $opt_p.1.12base.freq
if [ ! -s $opt_p.1.12base.freq ]; then
  echo "Error: build_lmer_table output error" >&2
  exit 100
fi

### 2. #create fasta file containing all kinds of repeats
### RepeatScout
RepeatScout -sequence $opt_i -freq $opt_p.1.12base.freq -output $opt_p.2.firstrepeats.lib
if [ ! -s $opt_p.2.firstrepeats.lib ]; then
  echo "Error: RepeatScout output error" >&2
  exit 100
fi

### 3. filter out short (<50bp) sequences
### RepeatScout: filter-stage-1.prl
cat $opt_p.2.firstrepeats.lib | filter-stage-1.prl > $opt_p.3.filtered.lib
if [ ! -s $opt_p.3.filtered.lib ]; then
  echo "Error: filter-stage-1.prl output error" >&2
  exit 100
fi


### 4. run RepeatMasker on your genome using filtered RepeatScout library
### RepeatMasker
if [ "$num_cpus" -gt "1" ]; then
  cmd_repeatmask="$cmd_repeatmask -pa $num_cpus"
fi
$cmd_repeatmask -lib $opt_p.3.filtered.lib -dir . $opt_i
#-nolow 
#output .cat .masked .out .ori.out tbl
if [ -s $fasta_name.out ]; then
  mv $fasta_name.out $opt_p.4.rm1.out
else
  echo "Error: RepeatMasker output error1" >&2
  exit 1
fi
if [ -s $fasta_name.cat ]; then
  mv $fasta_name.cat $opt_p.4.rm1.cat
fi
if [ -s $fasta_name.masked ]; then
  mv $fasta_name.masked $opt_p.4.rm1.masked
fi
if [ -s $fasta_name.ori.out ]; then
  mv $fasta_name.ori.out $opt_p.4.rm1.ori.out
fi
if [ -s $fasta_name.tbl ]; then
  mv $fasta_name.tbl $opt_p.4.rm1.tbl
fi
#The output file genome.fa.out will contain a column with the count of each repeat type.

### 5. filtering putative repeats by copy number
### RepeatScout: filter-stage-2.prl
#(--thres, by default, 10)
cat $opt_p.3.filtered.lib | filter-stage-2.prl --cat $opt_p.4.rm1.out --thresh $thres > $opt_p.5.secondfilter.lib
if [ ! -s $opt_p.5.secondfilter.lib ]; then
  echo "Error: filter-stage-2.prl output error" >&2
  exit 100
fi

### final RepeatMasker
if [ $run_final_repeatmasker -eq 1 ]; then
	echo "Info: run 2nd repeatMasker"
	$cmd_repeatmask -lib $opt_p.5.secondfilter.lib -dir . $opt_i
	if [ -s $fasta_name.out ]; then
	  mv $fasta_name.out $opt_p.6.rm2.out
	else
	  echo "Error: RepeatMasker output error2" >&2
	  exit 1
	fi
	if [ -s $fasta_name.cat ]; then
	  mv $fasta_name.cat $opt_p.6.rm2.cat
	fi
	if [ -s $fasta_name.masked ]; then
	  mv $fasta_name.masked $opt_p.6.rm2.masked
	fi
	if [ -s $fasta_name.ori.out ]; then
	  mv $fasta_name.ori.out $opt_p.6.rm2.ori.out
	fi
	if [ -s $fasta_name.tbl ]; then
	  mv $fasta_name.tbl $opt_p.6.rm2.tbl
	fi
fi	
#To view the length and frequency of repeats, open *.tbl.  Total length, %GC, and percent bases masked (i.e. percent in repeats) are shown. Also shown are number and total length occupied by different kinds of repeats. 


echo "Info: RepeatScout-pipe.sh done"
exit 0
