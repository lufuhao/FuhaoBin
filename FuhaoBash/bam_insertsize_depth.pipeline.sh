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

$0 --- Test library insert size and depth

Version: 20170214

Requirements:
    Linux: perl, md5sum, tail, wc, cut 
    Program: fastqc, flash, cutadapt, trimmomatic, bwa, samtools
    Script: fastq_checkid.pl, pcrdup, bamutils.pcrdup.insertsize.cluster.pl
            SizeCollectBin_luf.pl

Descriptions:
  Used to Test library insert size and depth
    1. FLASH to remove overlappable reads
    2. CutAdapt to remove adaptors
    3. Trimmomatic to remove low quality reads, Note phred33 only
    4. Test read file paired
    5. BWA aln to align each reads to reference
    6. Detect insert size, depth and redundancy

Options:
  -h    Print this help message
  -1    Read1.fq.gz
  -2    read2.fa.gz
  -r    Reference.fasta
  -x    BWA index
  -t    Number of threads, default: 1
  -a    Start step 1-6, default: 1
  -b    Stop step 1-6, default: 6
  -s    Comma delimited steps: 1,2,3,4,5,6
  -p    outprefix
  -z    Max insert size
  -bam  BAM file to run Step6 only
  -fo   FLASH options: -fo "-m Min_overlap -M Max_overlap -r read_length -f insert -s insert_stdev"
  -co   CutAdapt options
  -mo   Trimmomatic options: "ILLUMINACLIP:/path/to/seqprimer.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:150"

Example:
  $0 -1 R1.fq.gz -2 R2.fq.gz -r chr3b.fa -t 1 -z 6000 \
     -co "--times=3 --front=GATCTCTACCAGG --front=CCTGGTAGAG --overlap=8"
     -fo "-m 20 -M 251 -r 250 -f 800 -s 200"
     -mo "ILLUMINACLIP:$HOME/Hiseq.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:150"
     -p MyOutput.Miseq.test \
     -s 1,2,3,4,5,6

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



#################### Initializing ###################################
opt_r1=''
opt_r2=''
referenceseq=''
BwaIndex=''
numthreads=1
opt_bam=''
StepStart=1
StepEnd=6
RunDir=$PWD
flashoptions=''
outprefix='ReadOut'
customsteps=''
cutadaptoptions=''
trimoptions=''
maxinsert=0

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) opt_r1=$2;shift 2;;
    -2) opt_r2=$2;shift 2;;
    -r) referenceseq=$2; shift 2;;
    -fo) flashoptions=$2; shift 2;;
    -co) cutadaptoptions=$2; shift 2;;
    -mo) trimoptions=$2; shift 2;;
    -x) BwaIndex=$2; shift 2;;
    -bam) opt_bam=$2; shift 2;;
    -t) numthreads=$2;shift 2;;
    -s) customsteps=$2;shift 2;;
    -p) outprefix=$2;shift 2;;
    -a) StepStart=$2;shift 2;;
    -b) StepEnd=$2;shift 2;;
    -z) maxinsert=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" >&2; exit 100;;
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

File2Fullpath () {
	local FFinput=$1
	local FFsubinfo="SUF(File2Fullpath)"
	
	if [[ "$FFinput" =~ .*\/.* ]]; then
		if [[ "$FFinput" =~ ^\/ ]]; then
			echo $FFinput
		else
			echo $(cd $(dirname "$FFinput"); pwd)/$(basename "$FFinput")
		fi
	elif [[ "$FFinput" =~ ^[-0-9a-zA-Z\._]+$ ]]; then
		echo $(pwd)/$(basename "$FFinput")
	else
		echo ''
	fi
}
### Usage: RunFastQC $input_fastq $options
### Options: " -o /path/to/outputdir --nogroup"
### Global: $numthreads
### Dependancy: fastqc 
RunFastQC () {
	local RFfastainput=$1
	local RFfastoption=$2
	local RFsubinfo='SUF(RunFastQC)'
	
	fastqc $RFfastainput -f fastq -t $numthreads --noextract $RFfastoption
	ExitCode=$?
	if [ $ExitCode -ne 0 ]; then
		echo "${RFsubinfo}Error: fastqc running failed" >&2
		echo "CMD used fastqc $RFfastainput -f fastq -t $numthreads --noextract $RFfastoption" >&2
		echo 1;
	fi
	echo 0;
}
### BWA aln align fastq, sort and index bam
### Usage: RunFastQC $input_fastq $options
### Options: " -o /path/to/outputdir --nogroup"
### Global: $numthreads, $maxinsert
### Dependancy: fastqc 
RunBwaAln () {
	local RBAindex=$1
	local RBAfastqR1=$2
	local RBAfastqR2=$3
	local RBAoutbam=$4
	local RBAsubinfo="SUF(RunBwaAln)"
	
	bwa aln -t $numthreads $RBAindex $RBAfastqR1 > BWA_aln.temp.R1.sai
	if [ $? -ne 0 ] || [ ! -s "BWA_aln.temp.R1.sai" ]; then
		echo "${RBAsubinfo}Error: BWA aln R1 failed" >&2
		echo "CMD used: bwa aln -t $numthreads $RBAindex $RBAfastqR1 > BWA_aln.temp.R1.sai"
		return 1
	fi
	bwa aln -t $numthreads $RBAindex $RBAfastqR2 > BWA_aln.temp.R2.sai
	if [ $? -ne 0 ] || [ ! -s "BWA_aln.temp.R2.sai" ]; then
		echo "${RBAsubinfo}Error: BWA aln R2 failed" >&2
		echo "CMD used: bwa aln -t $numthreads $RBAindex $RBAfastqR2 > BWA_aln.temp.R2.sai"
		return 1
	fi
	bwa sampe -a $maxinsert $RBAindex "BWA_aln.temp.R1.sai" "BWA_aln.temp.R2.sai" $RBAfastqR1 $RBAfastqR2 | samtools view -bS -F 4 - > BWA_aln.temp.R1R2.bam
	if [ $? -ne 0 ] || [ ! -s "BWA_aln.temp.R1R2.bam" ]; then
		echo "${RBAsubinfo}Error: BWA sampe failed" >&2
		echo "CMD used: bwa sampe -a $maxinsert $RBAindex BWA_aln.temp.R1.sai BWA_aln.temp.R2.sai $RBAfastqR1 $RBAfastqR2 | samtools view -bS -F 4 - > BWA_aln.temp.R1R2.bam"
		return 1
	else
		rm "BWA_aln.temp.R1.sai" "BWA_aln.temp.R2.sai" > /dev/null 2>&1
	fi
	samtools sort BWA_aln.temp.R1R2.bam BWA_aln.temp.R1R2.st
	if [ $? -ne 0 ] || [ ! -s "BWA_aln.temp.R1R2.st.bam" ]; then
		echo "${RBAsubinfo}Error: samtools sort failed" >&2
		echo "CMD used: samtools sort -f BWA_aln.temp.R1R2.bam BWA_aln.temp.R1R2.st"
		return 1
	else
		rm "BWA_aln.temp.R1R2.bam" > /dev/null 2>&1
	fi
	
	mv "BWA_aln.temp.R1R2.st.bam" "$RBAoutbam"
	if [ $? -ne 0 ] || [ ! -s "$RBAoutbam" ]; then
		echo "${RBAsubinfo}Error: rename failed" >&2
		echo "CMD used: mv BWA_aln.temp.R1R2.st.bam $RBAoutbam"
		return 1
	fi
	samtools index $RBAoutbam
	if [ $? -ne 0 ]; then
		echo "${RBAsubinfo}Error: samtools index failed" >&2
		echo "CMD used: samtools index $RBAoutbam"
		return 1
	fi
	echo "${RBAsubinfo}### SAMtools FLAGSTAT ###"

	samtools flagstat $RBAoutbam
	
	return 0
}


### Decide which Step(s) to run
test_step1_flash=0;
test_step2_cutadapt=0;
test_step3_trimmomatic=0;
test_step4_pairness=0;
test_step5_bwa=0;
test_step6_insert=0;
if [ ! -z "$customsteps" ]; then
	allsteps=($(echo "$customsteps" | perl -ne 'chomp;@arr=split(/,/); foreach (@arr) {print $_, "\n" if (/^\d+$/);}'))
	for indstep in "${allsteps[@]}"; do
		case $indstep in
			1) test_step1_flash=1;;
			2) test_step2_cutadapt=1;;
			3) test_step3_trimmomatic=1;;
			4) test_step4_pairness=1;;
			5) test_step5_bwa=1;;
			6) test_step6_insert=1;;
			*) echo "Error: invalid -s $indstep. Should be a INT number between 0-6" >&2; exit 100;;
		esac
	done
else
	if [ $StepStart -le 1 ] && [ $StepEnd -ge 1 ]; then
		test_step1_flash=1;
	fi
	if [ $StepStart -le 2 ] && [ $StepEnd -ge 2 ]; then
		test_step2_cutadapt=1;
	fi
	if [ $StepStart -le 3 ] && [ $StepEnd -ge 3 ]; then
		test_step3_trimmomatic=1;
	fi
	if [ $StepStart -le 4 ] && [ $StepEnd -ge 4 ]; then
		test_step4_pairness=1;
	fi
	if [ $StepStart -le 5 ] && [ $StepEnd -ge 5 ]; then
		test_step5_bwa=1;
	fi
	if [ $StepStart -le 6 ] && [ $StepEnd -ge 6 ]; then
		test_step6_insert=1;
	fi
fi

if [ $test_step1_flash -eq 1 ] || [ $test_step2_cutadapt -eq 1 ] || [ $test_step3_trimmomatic -eq 1 ] || [ $test_step4_pairness -eq 1 ] || [ $test_step5_bwa -eq 1 ]; then
	if [ -z "$opt_r1" ] || [ ! -s "$opt_r1" ]; then
		echo "Error: invalid read R1 file" >&2
		exit 100
	fi
	if [ -z "$opt_r2" ] || [ ! -s "$opt_r2" ]; then
		echo "Error: invalid read R2 file" >&2
		exit 100
	fi
	opt_r1=$(File2Fullpath "$opt_r1")
	opt_r2=$(File2Fullpath "$opt_r2")
fi
if [ ! -z "$opt_bam" ] && [ -s "$opt_bam" ]; then
	opt_bam=$(File2Fullpath "$opt_bam")
	echo "Info: BAM input detected: $opt_bam"
	echo "Info: only running the final step"
	test_step1_flash=0;
	test_step2_cutadapt=0;
	test_step3_trimmomatic=0;
	test_step4_pairness=0;
	test_step5_bwa=0;
	test_step6_insert=1;
fi



#################### Command test ###################################

if [ ! -d "$RunDir/0.fastqc" ]; then
	if [ $test_step1_flash -eq 1 ] || [ $test_step2_cutadapt -eq 1 ] || [ $test_step3_trimmomatic -eq 1 ]; then
		if [ $(CmdExists 'fastqc') -ne 0 ]; then
			echo "Error: CMD 'fastqc' in PROGRAM 'FastQC' is required but not found.  Aborting..." >&2 
			exit 127
		fi
	fi
fi
if [ $test_step1_flash -eq 1 ]; then
	if [ $(CmdExists 'flash') -ne 0 ]; then
		echo "Error: CMD 'flash' in PROGRAM 'FLASH' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ -z "$flashoptions" ]; then
		echo "Error: FLASH options -fo need to be specified" >&2
		exit 100;
	fi
fi
if [ $test_step2_cutadapt -eq 1 ]; then
	if [ $(CmdExists 'cutadapt') -ne 0 ]; then
		echo "Error: CMD 'cutadapt' in PROGRAM 'CutAdapt' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ -z "$cutadaptoptions" ]; then
		echo "Error: cutadapt options -fo need to be specified" >&2
		exit 100;
	fi
fi
if [ $test_step3_trimmomatic -eq 1 ]; then
	if [ $(CmdExists 'trimmomatic') -ne 0 ]; then
		echo "Error: Script 'trimmomatic' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ -z "$trimoptions" ]; then
		echo "Error: cutadapt options -fo need to be specified" >&2
		exit 100;
	fi
fi
if [ $test_step4_pairness -eq 1 ]; then
	if [ $(CmdExists 'fastq_checkid.pl') -ne 0 ]; then
		echo "Error: Script 'fastq_checkid.pl' is required but not found.  Aborting..." >&2 
		exit 127
	fi
fi
if [ $test_step5_bwa -eq 1 ]; then
	if [ $(CmdExists 'bwa') -ne 0 ]; then
		echo "Error: CMD 'bwa' in PROGRAM 'BWA' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ ! -z "$BwaIndex" ]; then
		BwaIndex=$(File2Fullpath "$BwaIndex")
		echo "        Using existing BWA index $BwaIndex"
	fi
	if [ $maxinsert -gt 0 ] && [[ "$maxinsert" =~ ^[0-9]+$ ]]; then
		echo "" > /dev/null
	else
		echo "Error: -z maxinsert needed" >&2
		exit 100
	fi
	if [ $(CmdExists 'samtools') -ne 0 ]; then
		echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
		exit 127
	fi
fi
if [ $test_step6_insert -eq 1 ]; then
	if [ $(CmdExists 'pcrdup') -ne 0 ]; then
		echo "Error: Script 'pcrdup' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ $(CmdExists 'bamutils.pcrdup.insertsize.cluster.pl') -ne 0 ]; then
		echo "Error: Script 'bamutils.pcrdup.insertsize.cluster.pl' is required but not found.  Aborting..." >&2 
		exit 127
	fi
	if [ $(CmdExists 'SizeCollectBin_luf.pl') -ne 0 ]; then
		echo "Error: Script 'SizeCollectBin_luf.pl' is required but not found.  Aborting..." >&2 
		exit 127
	fi
fi


#################### Defaults #######################################
#outprefix=$(File2Fullpath $outprefix)



#################### Input and Output ###############################




### print input summary
echo "### Number threads: $numthreads"
if [ ! -z "$opt_r1" ] && [ -s "$opt_r1" ] && [ ! -z "$opt_r2" ] && [ -s "$opt_r2" ]; then
	echo "### Input FastQ R1: $opt_r1"
	echo "### Input FastQ R2: $opt_r2"
elif [ ! -z "$opt_bam" ] && [ -s "$opt_bam" ]; then
	echo "### Input BAM: $opt_bam"
else
	echo "Error: invalid input" >&2
	exit 100;
fi

echo "### Steps to run as follows:"
if [ $test_step1_flash -eq 1 ]; then
	echo "    Step1: FLASH: with option: $flashoptions"
fi
if [ $test_step2_cutadapt -eq 1 ]; then
	echo "    Step2: CutAdapt with option: $cutadaptoptions"
fi
if [ $test_step3_trimmomatic -eq 1 ]; then
	echo "    Step3: Trimmomatic with option $trimoptions"
fi
if [ $test_step4_pairness -eq 1 ]; then
	echo "    Step4: pairness"
fi
if [ $test_step5_bwa -eq 1 ]; then
	echo "    Step5: BWA with maxinsert: $maxinsert"
fi
if [ $test_step6_insert -eq 1 ]; then
	echo "    Step6: Analysis"
fi 






#################### Main ###########################################
if [ ! -z "$opt_r1" ] && [ -s "$opt_r1" ] && [ ! -z "$opt_r2" ] && [ -s "$opt_r2" ]; then
	echo "###### Step0: FastQC #####"
	echo "###### Step0: FastQC #####" >&2
	echo "Secret: to disable the primary FastQC step, mkdir $RunDir/0.fastqc"
	if [ ! -d $RunDir/0.fastqc ]; then
		mkdir -p $RunDir/0.fastqc
		cd $RunDir/0.fastqc
		if [ ! -d $RunDir/0.fastqc/group ]; then
			mkdir -p $RunDir/0.fastqc/group
			if [ $(RunFastQC "$opt_r1" "-o $RunDir/0.fastqc/group") -ne 0 ]; then
				echo "Error: fastqc R1 [grouped] failed: $opt_r1" >&2
				exit 100
			fi
			if [ $(RunFastQC "$opt_r2" "-o $RunDir/0.fastqc/group") -ne 0 ]; then
				echo "Error: fastqc R2 [grouped] failed: $opt_r2" >&2
				exit 100
			fi
		fi
		if [ ! -d $RunDir/0.fastqc/nogroup ]; then
			mkdir -p $RunDir/0.fastqc/nogroup
			if [ $(RunFastQC "$opt_r1" "-o $RunDir/0.fastqc/group --nogroup") -ne 0 ]; then
				echo "Error: fastqc R1 [no-grouped] failed: $opt_r1" >&2
				exit 100
			fi
			if [ $(RunFastQC "$opt_r2" "-o $RunDir/0.fastqc/group --nogroup") -ne 0 ]; then
				echo "Error: fastqc R2 [no-grouped] failed: $opt_r2" >&2
				exit 100
			fi
		fi
	fi
fi




if [ $test_step1_flash -eq 1 ]; then
	echo "###### Step1: FLASH #####"
	echo "###### Step1: FLASH #####" >&2
	CurrentPATH=$RunDir/1.flash
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	cd $CurrentPATH
	echo "Input: R1: $opt_r1"
	echo "Input: R1: $opt_r1" >&2
	echo "       R2: $opt_r2"
	echo "       R2: $opt_r2" >&2
# -m, --min-overlap=NUM
# -M, --max-overlap=NUM
# -r, --read-len=LEN
# -f, --fragment-len=LEN
# -s, --fragment-len-stddev=LEN
	flash $flashoptions --output-prefix=$outprefix -O -z -t $numthreads "$opt_r1" "$opt_r2"
	ExitCode=$?
	if [ $ExitCode -ne 0 ]; then
		echo "Error: flash running error" >&2
		echo "CMD used: flash $flashoptions --output-prefix=$outprefix -O -z -t $numthreads $opt_r1 $opt_r2" >&2
		exit $ExitCode;
	fi
	
	

	if [ -s "$outprefix.extendedFrags.fastq.gz" ]; then
		if [ $(RunFastQC "$outprefix.extendedFrags.fastq.gz" "--nogroup") -ne 0 ]; then
			echo "Error: fastqc merged paired [no-grouped] failed: $opt_r1" >&2
			exit 100
		fi
		zcat $outprefix.extendedFrags.fastq.gz | perl -ne 'chomp; if ($_=~/^\@(\S+)\s*/) {$ID=$1;}else{die "IDmatcherror\n";} $seq=<>; chomp $seq; $length=length($seq); <>; <>; print $ID."\t".$length."\n";' > $outprefix.extendedFrags.fastq.length
		if [ $? -ne 0 ]; then
			echo "Warnings: collect length failed: $outprefix.extendedFrags.fastq.gz" >&2
		fi
		if [ -s "$outprefix.extendedFrags.fastq.length" ]; then
			overlapedseqnum=$(wc -l < $outprefix.extendedFrags.fastq.length)
			echo "Info: Number of overlaped reads: $overlapedseqnum"
		
			SizeCollectBin_luf.pl $outprefix.extendedFrags.fastq.length 20 > $outprefix.extendedFrags.fasta.SizeBin
			if [ $? -ne 0 ] || [ ! -s "$outprefix.extendedFrags.fasta.SizeBin" ]; then
				echo "Warning: collect size bin failed" >&2
				echo "CMD used: SizeCollectBin_luf.pl $outprefix.extendedFrags.fastq.length 20 > $outprefix.extendedFrags.fasta.SizeBin" >&2
			fi
		
			### [OPTIONAL]: check high-occurance overlaped reads
			#fastq2fasta $outprefix.extendedFrags.fastq.gz > $outprefix.extendedFrags.fasta
			#cd-hit-est -i $outprefix.extendedFrags.fasta -o $outprefix.extendedFrags.cdhit -c 0.90 -n 8 -T 0 -r 1 -d 0 -M 30000
		
			### [OPTIONAL]: trimm overlap
			### Only useful when necessary, for instance. large proportion of reads are overlaped
			### Need: seqtk fastq_trim_overlap.pl
			#cat $outprefix.extendedFrags.fastq.length | cut -f 1 > $outprefix.extendedFrags.fastq.ID
			#seqtk subseq "$opt_r1" $outprefix.extendedFrags.fastq.ID > $outprefix.seqtk.R1.fq
			#seqtk subseq "$opt_r2" $outprefix.extendedFrags.fastq.ID > $outprefix.seqtk.R2.fq
			#fastq_trim_overlap.pl --flashfile $outprefix.extendedFrags.fastq.gz --forward $ID1.seqtk.R1.fq --reverse $outprefix.seqtk.R2.fq --outR1 $outprefix.nomerge.R1.fq --outR2 $outprefix.nomerge.R2.fq
		else
			echo "Warnings: collect length output error: $outprefix.extendedFrags.fastq.length" >&2
		fi
	fi
	
	opt_r1=$(File2Fullpath $outprefix.notCombined_1.fastq.gz)
	opt_r2=$(File2Fullpath $outprefix.notCombined_2.fastq.gz)
	if [ -s "$opt_r1" ] && [ -s "$opt_r2" ]; then
		if [ $(RunFastQC "$opt_r1" "--nogroup") -ne 0 ]; then
			echo "Error: fastqc R1 [no-grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(RunFastQC "$opt_r2" "--nogroup") -ne 0 ]; then
			echo "Error: fastqc R1 [no-grouped] failed: $opt_r2" >&2
			exit 100
		fi
	else
		echo "Error: flash output error" >&2
		echo "CMD used: flash $flashoptions --output-prefix=$outprefix -O -z -t $numthreads $opt_r1 $opt_r2" >&2
		exit 100;
	fi
	echo "Output: R1: $opt_r1"
	echo "Output: R1: $opt_r1" >&2
	echo "        R2: $opt_r2"
	echo "        R2: $opt_r2" >&2
fi



if [ $test_step2_cutadapt -eq 1 ]; then
	echo "###### Step2: CutAdapt #####"
	echo "###### Step2: CutAdapt #####" >&2
	CurrentPATH=$RunDir/2.cutadapt
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	cd $CurrentPATH
	echo "Input: R1: $opt_r1"
	echo "Input: R1: $opt_r1" >&2
	echo "       R2: $opt_r2"
	echo "       R2: $opt_r2" >&2
	cutadapt --format=fastq $cutadaptoptions -o $outprefix.cutadapt_R1.fq.gz "$opt_r1"
	if [ $? -ne 0 ] || [ ! -s "$outprefix.cutadapt_R1.fq.gz" ]; then
		echo "Error: cutadapt R1 runnng failed" >&2
		echo "CMD used cutadapt --format=fastq $cutadaptoptions -o $outprefix.cutadapt_R1.fq.gz $opt_r1" >&2
		exit 100
	fi
	cutadapt --format=fastq $cutadaptoptions -o $outprefix.cutadapt_R2.fq.gz "$opt_r2"
	if [ $? -ne 0 ] || [ ! -s "$outprefix.cutadapt_R2.fq.gz" ]; then
		echo "Error: cutadapt R2 runnng failed" >&2
		echo "CMD used cutadapt --format=fastq $cutadaptoptions -o $outprefix.cutadapt_R2.fq.gz $opt_r2" >&2
		exit 100
	fi
	
	opt_r1=$(File2Fullpath "$outprefix.cutadapt_R1.fq.gz")
	opt_r2=$(File2Fullpath "$outprefix.cutadapt_R2.fq.gz")
	if [ -s "$opt_r1" ] && [ -s "$opt_r1" ]; then
		if [ $(RunFastQC "$opt_r1" "--nogroup") -ne 0 ]; then
			echo "Error: fastqc R1 [no-grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(RunFastQC "$opt_r2" "--nogroup") -ne 0 ]; then
			echo "Error: fastqc R2 [no-grouped] failed: $opt_r2" >&2
			exit 100
		fi
	fi
	echo "Output: R1: $opt_r1"
	echo "Output: R1: $opt_r1" >&2
	echo "        R2: $opt_r2"
	echo "        R2: $opt_r2" >&2
fi



if [ $test_step3_trimmomatic -eq 1 ]; then
	echo "###### Step3: trimmomatic #####"
	echo "###### Step3: trimmomatic #####" >&2
	CurrentPATH=$RunDir/3.trimmomatic
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	mkdir $CurrentPATH/nogroup
	mkdir $CurrentPATH/group
	cd $CurrentPATH
	echo "Input: R1: $opt_r1"
	echo "Input: R1: $opt_r1" >&2
	echo "       R2: $opt_r2"
	echo "       R2: $opt_r2" >&2

	trimmomatic PE -threads $numthreads -phred33 -trimlog $outprefix.R1_R2.trimlog $opt_r1 $opt_r2 $outprefix.R1.trim.fq.gz $outprefix.R1.unpaired.fq.gz $outprefix.R2.trim.fq.gz $outprefix.R2.unpaired.fq.gz $trimoptions
	if [ $? -ne 0 ] || [ ! -s "$outprefix.R1.trim.fq.gz" ] || [ ! -s "$outprefix.R2.trim.fq.gz" ]; then
		echo "Error: trimmomatic error" >&2
		echo "CMD used: trimmomatic PE -threads $numthreads -phred33 -trimlog $outprefix.R1_R2.trimlog $opt_r1 $opt_r2 $outprefix.R1.trim.fq.gz $outprefix.R1.unpaired.fq.gz $outprefix.R2.trim.fq.gz $outprefix.R2.unpaired.fq.gz $trimoptions" >&2
		exit 100
	fi
	opt_r1=$(File2Fullpath "$outprefix.R1.trim.fq.gz")
	opt_r2=$(File2Fullpath "$outprefix.R2.trim.fq.gz")
	if [ -s "$opt_r1" ] && [ -s "$opt_r2" ]; then
		if [ $(RunFastQC "$opt_r1" "-o $CurrentPATH/nogroup --nogroup") -ne 0 ]; then
			echo "Error: fastqc R1 [no-grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(RunFastQC "$opt_r2" "-o $CurrentPATH/nogroup --nogroup") -ne 0 ]; then
			echo "Error: fastqc R1 [no-grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(RunFastQC "$opt_r1" "-o $CurrentPATH/group") -ne 0 ]; then
			echo "Error: fastqc R1 [grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(RunFastQC "$opt_r2" "-o $CurrentPATH/group") -ne 0 ]; then
			echo "Error: fastqc R1 [grouped] failed: $opt_r1" >&2
			exit 100
		fi
		if [ $(CmdExists 'md5sum') -eq 0 ]; then
			md5sum "$opt_r1" > "$opt_r1.md5"
			md5sum "$opt_r2" > "$opt_r1.md5"
		else
			echo "Warnings: md5sum command not found, skipping" >&2 
		fi
	else
		echo "Error: trimmomatic output error" >&2
		exit 100
	fi
	
	echo "Output: R1: $opt_r1"
	echo "Output: R1: $opt_r1" >&2
	echo "        R2: $opt_r2"
	echo "        R2: $opt_r2" >&2
fi



if [ $test_step4_pairness -eq 1 ]; then
	echo "###### Step4: pairness #####"
	echo "###### Step4: pairnees #####" >&2
	CurrentPATH=$RunDir/4.pairnees
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	cd $CurrentPATH
	echo "Input: R1: $opt_r1"
	echo "Input: R1: $opt_r1" >&2
	echo "       R2: $opt_r2"
	echo "       R2: $opt_r2" >&2
	fastq_checkid.pl "$opt_r1" "$opt_r2" '\@(\S+)\s*\S*'
	if [ $? -eq 0 ]; then
	  echo "Info: paired fastq: $opt_r1 $opt_r2"
	  touch ${$outprefix}_paired
	else
	  echo "Error: unpaired fastq: $opt_r1 $opt_r2"
	  touch ${$outprefix}_unpaired
	  exit 100
	fi
	
	echo "Output: R1: $opt_r1"
	echo "Output: R1: $opt_r1" >&2
	echo "        R2: $opt_r2"
	echo "        R2: $opt_r2" >&2
fi



if [ $test_step5_bwa -eq 1 ]; then
	echo "###### Step5: BWA #####"
	echo "###### Step5: BWA #####" >&2
	CurrentPATH=$RunDir/5.bwa
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	cd $CurrentPATH
	echo "Input: R1: $opt_r1"
	echo "Input: R1: $opt_r1" >&2
	echo "       R2: $opt_r2"
	echo "       R2: $opt_r2" >&2
	opt_bam=$outprefix.sort.bam
	out_bam=$(File2Fullpath "$opt_bam")
	
	if [ -z "$BwaIndex" ]; then
		bwa index -a is -p $outprefix.bwaindex $referenceseq
		if [ $? -ne 0 ]; then
			echo "Error: BWA index error" >&2
			echo "CMD used: bwa index -a is -p $outprefix.bwaindex $referenceseq" >&2
			exit 100
		fi
		BwaIndex=$(File2Fullpath "$outprefix.bwaindex")
	fi
	
	if RunBwaAln $BwaIndex $opt_r1 $opt_r2 $opt_bam; then
		echo "Info: BWA running is successful, BAM: $opt_bam" >/dev/null
	else
		echo "Error: BWA aln mapping failed" >&2
		exit 100
	fi
	
	echo "Output: BAM: $opt_bam"
	echo "Output: BAM: $opt_bam" >&2
fi



if [ $test_step6_insert -eq 1 ]; then
	echo "###### Step6: Analysis #####"
	echo "###### Step6: Analysis #####" >&2
	if [ -z "$opt_bam" ] || [ ! -s "$opt_bam" ]; then
		echo "Error: BAM file not found" >&2
		exit 100
	fi
	CurrentPATH=$RunDir/6.insertsize
	if [ -d $CurrentPATH ]; then
		rm -rf "$CurrentPATH" > /dev/null 2>&1
	fi
	mkdir -p $CurrentPATH
	cd $CurrentPATH
	echo "Input: BAM: $opt_bam"
	echo "Input: BAM: $opt_bam" >&2
	pwd
	pcrdup -counts $outprefix.sort.pcrdup $opt_bam
	if [ $? -ne 0 ] || [ ! -s "$outprefix.sort.pcrdup" ]; then
		echo "Error: pcrdup failed" >&2
		exit 100
	fi
	bamutils.pcrdup.insertsize.cluster.pl $outprefix.sort.pcrdup 20 20 > $outprefix.sort.pcrdup.rmdup
	if [ $? -ne 0 ] || [ ! -s "$outprefix.sort.pcrdup.rmdup" ]; then
		echo "Error: pcrdup insertsize cluster failed" >&2
		exit 100
	fi
	echo "Total unique"
	tail -n +6 $outprefix.sort.pcrdup.rmdup | wc -l
	tail -n +6 $outprefix.sort.pcrdup.rmdup | perl -lane 'print "$F[2]\t$F[3]";' > $outprefix.sort.pcrdup.rmdup.col23
	tail -n +6 $outprefix.sort.pcrdup.rmdup | cut -f 4 > $outprefix.sort.pcrdup.rmdup.depth
	tail -n +6 $outprefix.sort.pcrdup.rmdup | cut -f 2 > $outprefix.sort.pcrdup.rmdup.pos
	cat $outprefix.sort.pcrdup.rmdup.depth | perl -lne '$sum+=$_;$count++;END {print "AverageDepth: ". $sum/$count. " ";}'
	tail -n +6  $outprefix.sort.pcrdup.rmdup | cut -f 3 > $outprefix.sort.pcrdup.rmdup.insert

	echo "pos 10M: $outprefix.sort.pcrdup.rmdup.pos.10Msizebin"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.pos 10000000 > $outprefix.sort.pcrdup.rmdup.pos.10Msizebin
	echo "pos 100K: $outprefix.sort.pcrdup.rmdup.pos.100Ksizebin"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.pos 100000 > $outprefix.sort.pcrdup.rmdup.pos.100Ksizebin

	echo "Unique per pos: $outprefix.sort.pcrdup.rmdup.insert.bin2000"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.insert 2000 > $outprefix.sort.pcrdup.rmdup.insert.bin2000
	echo "All perl pos: $outprefix.sort.pcrdup.rmdup.allpos.bin2000"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.col23 2000 > $outprefix.sort.pcrdup.rmdup.allpos.bin2000
	echo "duplicates : $outprefix.sort.pcrdup.rmdup.depth.bin2"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.depth 2 > $outprefix.sort.pcrdup.rmdup.depth.bin2
	echo "duplicates : $outprefix.sort.pcrdup.rmdup.depth.bin5"
	SizeCollectBin_luf.pl $outprefix.sort.pcrdup.rmdup.depth 5 > $outprefix.sort.pcrdup.rmdup.depth.bin5
fi


echo -e "\n######################\nProgram $ProgramName END\n######################\n"
exit 0;
