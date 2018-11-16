#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

run_splitNreads.sh

Version: 20150508

Require:
  samtools
  bam_splitNreads.pl

Descriptions:
  split cigar N into complementart parts using bam_splitNreads.pl
    bam_splitNreads.pl
    samtools calmd
    samtools sort
    samtools index

Options:
  -h    Print this help message

Example:
  run_splitNreads.sh 1.bam_in 2.reference 3.bam_out

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
[ -z "$2" ] && help
[ -z "$3" ] && help
[ "$1" = "-h" ] && help
[ "$1" = "--help" ] && help
#################### Defaults #######################################
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH



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
if [ $(CmdExists 'bam_splitNreads.pl') -eq 1 ]; then
	echo "Error: script 'bam_splitNreads.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'samtools') -eq 1 ]; then
	echo "Error: CMD 'samtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



###Input and output
input=$1
reference=$2
output=$3



###
input_name=${input##*/}
input_base=${input_name%.*}
###Step1 splitN
bam_splitNreads.pl $input | samtools view -b -S -h - > $input_base.splitN.bam
if [ $? -eq 0 ] && [ -s $input_base.splitN.bam ]; then
	echo "Step1 splitN: pass"
else
	echo "Step1 splitN: fail" > /dev/stderr
	exit 1
fi


###Step2 calmd
samtools calmd -b $input_base.splitN.bam $reference > $input_base.splitN.calmd.bam
if [ $? -eq 0 ] && [ -s $input_base.splitN.calmd.bam ]; then
	echo "Step2 calmd: pass"
else
	echo "Step2 calmd: fail" > /dev/stderr
	exit 1
fi



###Step3 sort
samtools sort -f $input_base.splitN.calmd.bam $output
if [ $? -eq 0 ] && [ -s $output ]; then
	echo "Step3 sort: pass"
else
	echo "Step3 sort: fail" > /dev/stderr
	exit 1
fi



### Step4 index
samtools index $output
if [ $? -eq 0 ] && [ -s $output.bai ]; then
	echo "Step4 index: pass"
#	rm $input_base.splitN.bam $input_base.splitN.calmd.bam
else
	echo "Step4 index: fail" > /dev/stderr
	exit 1
fi

###End
exit 0
