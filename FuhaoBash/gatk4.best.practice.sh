#!/bin/bash
### Exit if command fails
#set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
#set -o nounset
### Script Root
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
### MachType
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

###echo color
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
#RED='\033[0;31m'
#NC='\033[0m' # No Color
#printf "I ${RED}love${NC} Stack Overflow\n"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: v20200320

Requirements:
	GATK
	tabix, bgzip
	bcftools

Descriptions:
    1: MarkDuplicate
    2: Set up SAM tag
    3： Add read group
    4: Collect known SNPs
    5: BaseRecalibrator
    6: ApplyBQSR
    7: HaplotypeCaller
    8: CombineGVCFs
    9： GenotypeGVCFs

Options:
  -h    Print this help message
  -f    Reference fasta in flat txt
  -b    BAM files (sorted and index), delimited by comma
  -c    [Optional] Known SNPs in VCF format
  -r    [Optional] Read Groups delimited by comma
  -p    Out prefix for each BAM, delimited by comma
  -d    Output directory [default: PWD]
  -t    Number of threads, default: 1
  -m    Max memory for JAVA -XmxNNg option [INT, default: 4]

Example:
  gatk4.best.practice.sh -f $PWD/A01_70.fa -p AA,DD -m 8 -r AA,DD \
      -b $PWD/AAsubgen.AAseq.sort.A01.bam,$PWD/AAsubgen.DDseq.sort.A01.bam
      
  ** Note: make sure all the reference fasta in BAM exist in reference
               or HaplotyprClaaer report error


Author:
  Fu-Hao Lu
  Professor, PhD
  State Key Labortory of Crop Stress Adaptation and Improvement
  College of Life Science
  Jinming Campus, Henan University
  Kaifeng 475004, P.R.China
  E-mail: lufuhao@henu.edu.cn
HELP
exit 2
}
[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_f=""
opt_t=1
opt_m=4
declare -a bamArr=();
declare -a outPfx=();
declare -a readgroup=();
opt_out_dir=$PWD

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -f) opt_f=$2;shift 2;;
    -b) bamArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -p) outPfx=($(echo $2 | tr ',' "\n"));shift 2;;
    -r) readgroup=($(echo $2 | tr ',' "\n"));shift 2;;
    -c) vcfArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -d) opt_out_dir=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    return 1
  fi
}
printMsg () {
	echo "$1"
	echo "$1" >&2
}




#################### Command test ###################################
CmdExists 'bgzip'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'bgzip' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'tabix'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'tabix' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bcftools'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'bcftools' in PROGRAM 'BCFtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'gatk'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'gatk' in PROGRAM 'GATK' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
step=0;
opt_jopt=" -Xmx${opt_m}G -XX:ParallelGCThreads=$opt_t "
vcfOptions=""
if [ ${#vcfArr[@]} -gt 0 ]; then
	for indVcf in ${vcfArr[@]}; do
		vcfOptions=" $vcfOptions --known-sites $indVcf "
	done
fi



#################### Input and Output ###############################
if [ -z "${opt_f}" ] || [ ! -s "$opt_f" ]; then
	echo "Error: invalid reference fasta file" >&2
	exit 100
fi
if [ ! -s "${opt_f%.*}.dict" ]; then
	gatk CreateSequenceDictionary -R $opt_f -O ${opt_f%.*}.dict
	if [ $? -ne 0 ] || [ ! -s "${opt_f%.*}.dict" ]; then
		echo "Error: gatk CreateSequenceDictionary: $opt_f" >&2
		exit 100
	else
		printMsg "    Dict: done $opt_f"
	fi
fi
if [ ${#bamArr[@]} -eq 0 ]; then
	echo "Error: empty BAM files" >&2
	exit 100
else
	for indbam in "${bamArr[@]}"; do
		if [ ! -s "$indbam" ]; then
			echo "Error: invalid BAM: $indbam" >&2
			exit 100
		fi
	done
	if [ ${#outPfx[@]} -ne ${#bamArr[@]} ]; then
		echo "Error: unequal number between BAM files and out prefix (-p)" >&2
		exit 100
	fi
fi
if [ ${#readgroup[@]} -gt 0 ]; then
	if [ ${#readgroup[@]} -ne ${#bamArr[@]} ]; then
		echo "Error: unequal number between BAM files and read group (-r)" >&2
		exit 100
	fi
fi



#################### Main ###########################################



### step1: MarkDuplicate
((step++))
printMsg "Step${step}: MarkDuplicate"
outDir="$opt_out_dir/${step}.MarkDuplicates"
declare -a bamArr2=()
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
for ((i=0;i<${#outPfx[@]}; i++)); do
	inBam=${bamArr[$i]}
	outBam="$outDir/${outPfx[$i]}.${step}.markdup.bam"
	if [ ! -s "$outBam" ]; then
		printMsg "    BAM: $inBam"
		printMsg "    pfx: ${outPfx[$i]}"
		printMsg "    out: $outBam"
		gatk --java-options "$opt_jopt" MarkDuplicates -I $inBam -M ${outPfx[$i]}.${step}.markdup_metrics.txt -O $outBam
		if [ $? -ne 0 ] || [ ! -s $outBam ]; then
			echo "Error: MarkDuplicates: $inBam" >&2
			exit 100
		else
			printMsg "    BAM: $inBam MarkDuplicates done"
		fi
	else
		printMsg "Warnings: using existing deduplicated BAM: $outBam"
	fi
	bamArr2+=( "$outBam" )
	if [ ! -s "$outBam.bai" ]; then
		samtools index $outBam
		if [ $? -ne 0 ] || [ ! -s "$outBam.bai" ]; then
			echo "Error: index: $outBam" >&2
			exit 100
		else
			printMsg "    BAM: $outBam index done"
		fi
	fi
done
bamArr=("${bamArr2[@]}"); bamArr2=()



### step2: set up SAM tag: NM, MD, UQ
((step++))
printMsg "Step${step}: tagSAM"
outDir="$opt_out_dir/${step}.tagSAM"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
for ((i=0;i<${#outPfx[@]}; i++)); do
	inBam=${bamArr[$i]}
	outBam="$outDir/${outPfx[$i]}.${step}.tag.bam"
	if [ ! -s "$outBam" ]; then
		printMsg "    BAM: $inBam"
		printMsg "    out: $outBam"
		gatk --java-options "$opt_jopt" SetNmMdAndUqTags -I $inBam -O $outBam --CREATE_INDEX true --CREATE_MD5_FILE true --REFERENCE_SEQUENCE $opt_f 
		if [ $? -ne 0 ] || [ ! -s $outBam ]; then
			echo "Error: tagSAM: $inBam" >&2
			exit 100
		else
			printMsg "    BAM: $inBam tagSAM done"
		fi
	else
		printMsg "Warnings: using existing deduplicated BAM: $outBam"
	fi
	bamArr2+=("$outBam")
done
bamArr=("${bamArr2[@]}"); bamArr2=()



### Step 3: collect known SNPs
if [ ${#readgroup[@]} -gt 0 ]; then
	((step++))
	printMsg "Step${step}: add read group"
	outDir="$opt_out_dir/${step}.readgroup"
	if [ ! -d $outDir ]; then
		mkdir -p $outDir
	fi
	cd $outDir
	for ((i=0;i<${#outPfx[@]}; i++)); do
		inBam=${bamArr[$i]}
		outBam="$outDir/${outPfx[$i]}.${step}.RG.bam"
		if [ ! -s "$outBam" ]; then
			printMsg "    BAM: $inBam"
			printMsg "    out: $outBam"
			gatk --java-options "$opt_jopt" AddOrReplaceReadGroups --INPUT $inBam --OUTPUT $outBam --CREATE_INDEX true --REFERENCE_SEQUENCE $opt_f --RGPU "${readgroup[$i]}" --RGID "${readgroup[$i]}" --RGLB "${readgroup[$i]}" --RGPL "Illumina" --RGSM "${readgroup[$i]}"
			if [ $? -ne 0 ] || [ ! -s $outBam ]; then
				echo "Error: tagSAM: $inBam" >&2
				exit 100
			else
				printMsg "    BAM: $inBam tagSAM done"
			fi
		else
			printMsg "Warnings: using existing deduplicated BAM: $outBam"
		fi
		bamArr2+=("$outBam")
	done
	bamArr=("${bamArr2[@]}"); bamArr2=()
fi



### Step 4: collect known SNPs
((step++))
printMsg "Step${step}: VCF"
outDir="$opt_out_dir/${step}.vcf"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
if [ ${#vcfArr[@]} -eq 0 ]; then
	rawvcf="$outDir/merge.raw.vcf"
	if [ -e "bam.list" ]; then
		rm -rf bam.list
	fi
	for ((i=0;i<${#bamArr[@]}; i++)); do
		echo "${bamArr[$i]}" >> bam.list
	done
	if [ ! -s "$rawvcf" ];then
		bcftools mpileup -b bam.list --fasta-ref $opt_f | bcftools call -mv -o $rawvcf
		if [ $? -ne 0 ] || [ ! -s $rawvcf ]; then
			echo "Error: bcftools mpileup: $inBam" >&2
			exit 100
		else
			printMsg "    BAM: $inBam BaseRecalibrator done"
		fi
		cat $rawvcf | bgzip > $rawvcf.gz
		tabix -pvcf $rawvcf.gz
	fi

#	vcfutils.pl varFilter -D9999 HS0674.var.vcf >HS0674.varFilter.vcf
#		-D INT    maximum read depth [10000000]
#	bgzip HS0674.varFilter.vcf
#	tabix -pvcf HS0674.varFilter.vcf.gz

	vcfOptions=" ${vcfOptions} --known-sites $rawvcf.gz "
fi


### Step 5: BaseRecalibrator
((step++))
printMsg "Step${step}: BaseRecalibrator"
outDir="$opt_out_dir/${step}.BaseRecalibrator"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
declare -a recalArr=()
for ((i=0;i<${#outPfx[@]}; i++)); do
	inBam=${bamArr[$i]}
	outRecal="$outDir/${outPfx[$i]}.${step}.recal.table"
	printMsg "    BAM: $inBam"
	printMsg "    out: $outRecal"
	if [ ! -s "$outRecal" ]; then
		time gatk --java-options "$opt_jopt" BaseRecalibrator -R $opt_f -I $inBam $vcfOptions -O $outRecal
		if [ $? -ne 0 ] || [ ! -s $outRecal ]; then
			echo "Error: BaseRecalibrator: $inBam" >&2
			exit 100
		else
			printMsg "    BAM: $inBam BaseRecalibrator done"
		fi
	fi
	recalArr+=("$outRecal")
done
if [ ${#recalArr[@]} -ne ${#bamArr[@]} ] || [ ${#recalArr[@]} -eq 0 ]; then
	echo "Error: BaseRecalibrator error " >&2
	exit 100
fi



### Step 6: ApplyBQSR
((step++))
printMsg "Step${step}: ApplyBQSR"
outDir="$opt_out_dir/${step}.ApplyBQSR"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
declare -a bamArr2=()
for ((i=0;i<${#bamArr[@]}; i++)); do
	inBam=${bamArr[$i]}
	outbam="$outDir/${outPfx[$i]}.${step}.bqsr.bam"
	printMsg "    BAM: $inBam"
	printMsg "    out: $outBam"
	if [ ! -s "$outbam" ]; then
		time gatk --java-options "$opt_jopt" ApplyBQSR --bqsr-recal-file ${recalArr[$i]} -R $opt_f -I $inBam -O $outbam
		if [ $? -ne 0 ] || [ ! -s $outbam ]; then
			echo "Error: ApplyBQSR: $inBam" >&2
			exit 100
		else
			printMsg "    BAM: $inBam ApplyBQSR done"
		fi
	fi
	bamArr2+=("$outbam")
done
bamArr=("${bamArr2[@]}");  bamArr2=()



### Step 7: HaplotypeCaller
((step++))
printMsg "Step${step}: HaplotypeCaller"
outDir="$opt_out_dir/${step}.HaplotypeCaller"
sample_gvcfs=""
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
for ((i=0;i<${#bamArr[@]}; i++)); do
	inBam=${bamArr[$i]}
	outVcf="$outDir/${outPfx[$i]}.${step}.HaplotypeCaller.gvcf.gz"
	outBam="$outDir/${outPfx[$i]}.${step}.HaplotypeCaller.bam"
	if [ ! -s "$outVcf" ]; then
		printMsg "    BAM: $inBam"
		printMsg "    out: $outBam"
		printMsg "    out: $outVcf"
#--emit-ref-confidence,-ERC <ReferenceConfidenceMode>
#                              Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature) 
#                              Default value: NONE. Possible values: {NONE, BP_RESOLUTION, GVCF} 
		time gatk --java-options "$opt_jopt" HaplotypeCaller --emit-ref-confidence GVCF -R $opt_f -I $inBam -O $outVcf -bamout $outBam
		if [ $? -ne 0 ] || [ ! -s $outVcf ]; then
			echo "Error: HaplotypeCaller: $inBam" >&2
			echo "EMD used: time gatk --java-options $opt_jopt HaplotypeCaller --emit-ref-confidence GVCF -R $opt_f -I $inBam -O $outVcf -bamout $outBam"
			exit 100
		else
			printMsg "    BAM: $inBam HaplotypeCaller done"
		fi
	fi
	sample_gvcfs=" $sample_gvcfs -V $outVcf "
done



### Step 8: CombineGVCFs
((step++))
printMsg "Step${step}: CombineGVCFs"
outDir="$opt_out_dir/${step}.CombineGVCFs"
merged_gvcfs="$outDir/merged.gvcf.gz"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
if [ ! -s "merged_gvcfs" ]; then
	time gatk --java-options "$opt_jopt" CombineGVCFs -R $opt_f ${sample_gvcfs} -O $merged_gvcfs
	if [ $? -ne 0 ] || [ ! -s $merged_gvcfs ]; then
		echo "Error: CombineGVFs" >&2
		exit 100
	else
		printMsg "    CombineGVFs done"
	fi
fi



### Step 9: GenotypeGVCFs
((step++))
printMsg "Step${step}: GenotypeGVCFs"
outDir="$opt_out_dir/${step}.GenotypeGVCFs"
merged_vcf="$outDir/merged.vcf.gz"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
if [ ! -s "$merged_vcf" ]; then
	time gatk --java-options "$opt_jopt" GenotypeGVCFs -R $opt_f -V $merged_gvcfs -O $merged_vcf
	if [ $? -ne 0 ] || [ ! -s $merged_vcf ]; then
		echo "Error: GenotypeGVCFs" >&2
		exit 100
	else
		printMsg "    GenotypeGVCFs done"
	fi
fi


:<<EOM
((step++))
printMsg "Step${step}: VariantRecalibrator"
outDir="$opt_out_dir/${step}.VariantRecalibrator"
merged_vcf="$outDir/merged.vcf.gz"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
gatk --java-options "$opt_jopt" VariantRecalibrator -R $opt_f -V $merged_vcf -mode SNP \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -O ${outDir}/sample.recal --tranches-file ${outDir}/sample.tranches --rscript-file ${outDir}/sample.plots.R
if [ $? -ne 0 ] || [ ! -s $merged_gvcfs ]; then
	printMsg "    VariantRecalibrator"
else
	echo "Error: VariantRecalibrator"
	exit 100
fi



gatk --java-options "$opt_jopt" ApplyVQSR \
    -R ${REF_FA} \
    -V ./data/sample.g.vcf.gz \
    -O ./data/sample.recalibrated.g.vcf.gz \
    -OVI \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ./data/sample.tranches \
    --recal-file ./data/sample.recal \
    -mode SNP
EOM


exit 0
