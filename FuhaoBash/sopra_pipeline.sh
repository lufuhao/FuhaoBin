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

Version: 20160730

Requirements:
	perl && File::Spec

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
if [ $(CmdExists '') -eq 1 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi

s_prep_contigAseq_v1.4.6.pl
bowtie-build
bowtie
#################### Defaults #######################################

contigs=
readfas=
readname=${readfas##*/}
readbase=${readname%.*}
outdir=
copraroot=
numlinks=5
$insertsize
#################### Input and Output ###############################




#################### Main ###########################################
if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	exit 1
fi

### Step1: generate the new files contigs_sopra.fasta and srr_sopra.fasta
###change into the sopra with prebuilt contigs directory
###cd source_codes_v1.4.6/SOPRA_with_prebuilt_contigs/
perl $copraroot/s_prep_contigAseq_v1.4.6.pl -contig $contigs -mate $readfas -a $outdir/
if [ $? ne 0 ] || [ ! -s $outdir/contigs_sopra.fasta ] || [ ! -s ${readbase}_sopra.fasta ]; then
	echo "Error: Prepare contigs error\n" >&2
	exit 1
fi

### Step2: manually align the reads srr_sopra.fasta against the contigs
#We do this using bowtie using the -v 0 -f --sam -m 1 options (no mismatches, fasta format sam output and only report one reportable alignment, suppress alignments where multiple alignments are valid. The latter is not necessary as SOPRA would do this as well.)
cd $outdir
bowtie-build $outdir/contigs_sopra.fasta SOPRA
bowtie -v 0 -m 1 -f --sam SOPRA $outdir/${readbase}_sopra.fasta  > $outdir/mysam_mate_illu1
if [ $? -ne 0 ] || [ ! -s $outdir/mysam_mate_illu1 ]; then
	echo "Error: bowtie mapping error" >&2
	exit 1
fi

### Step 3: Parse sam
###This needs around 20GB of RAM
perl $copraroot/s_parse_sam_v1.4.6.pl -sam mysam_mate_illu1 -a $outdir/
if [ $? -ne 0 ] || [ ! -s $outdir/mysam_mate_illu1_parsed ]; then
	echo "Error: bowtie mapping error" >&2
	exit 1
fi


### Step 4: 
#Now we just give the insert size (200) of the mate pairs
perl $copraroot/s_read_parsed_sam_v1.4.6.pl  -c 5 -parsed mysam_mate_illu1_parsed -d $insertsize -a $outdir/
if [ $? -ne 0 ] || [ ! -s $outdir/mysam_mate_illu1_parsed ]; then
	echo "Error: bowtie mapping error" >&2
	exit 1
fi

### Step 5: scaffold
perl $copraroot/s_scaf_v1.4.6.pl -w $numlinks -o orientdistinfo_c5 -a $outdir/

### Step 6: Evaluation
#running fac.pl (You can get it from this email http://www.bcgsc.ca/pipermail/abyss-users/2009-September/000330.html)
contigs_stats.pl -i scaffolds_h2.2_L150_w5.fasta




# STEP 3
perl s_read_parsed_sam_v1.4.5.pl -parsed ${OUTDIR}/readstocontigs.sam_parsed -d ${INSERTSIZE} -a ${OUTDIR} -c ${C} # -e ${E}

# STEP 4: scaffold assembly
perl s_scaf_v1.4.5.pl -o ${OUTDIR}/orientdistinfo_c${C} -a ${OUTDIR} # -w ${W} -h ${H}

cp ${OUTDIR}/scaffold* .
