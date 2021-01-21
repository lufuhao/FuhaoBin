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
opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) FastQR1Arr=($(echo $2 | tr  "\n"));shift 2;;
    -t) opt_t=$2;shift 2;;
    -1) seq_rfn=(${seq_rfn[@]} "$2");shift 2;;
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

if [[ $(CmdExists 'mum.stat') -eq 1 ]]; then
	echo "Error: script 'mum.stat' is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################




#################### Input and Output ###############################
opt_out_dir=""
bamArr=() #sorted and index
outPfx=()
readGroup=()
opt_t=1
opt_m=10
opt_f="xxx.fa"
vcfArr=()
step=0;
opt_jopt="-Xmx10G -XX:ParallelGCThreads=4"

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
	printMsg "    BAM: $inBam"
	gatk MarkDuplicates -I $inBam -M ${outPfx[$i]}.${step}.markdup_metrics.txt -O $outBam
	if [ $? -ne 0 ] || [ ! -s $outBam ]; then
		bamArr2+=( "$outBam" )
		printMsg "    BAM: $inBam MarkDuplicates done"
	else
		echo "Error: MarkDuplicates: $inBam"
		exit 100
	fi
	samtools index $outBam
	if [ $? -ne 0] || [ ! -s $outBam.bai ]; then
		printMsg "    BAM: $outBam index done"
	else
		echo "Error: index: $outBam" >&2
		exit 100
	fi
done
bamArr=("${bamArr2[@]}"); bamArr2=()


### step2: set up SAM tag: NM, MD, UQ
### NM?????? ??????????????????MD??????????UQ??Phred ??????????????????
((step++))
printMsg "Step${step}: tagSAM"
outDir="$opt_out_dir/${step}.tagSAM"
bamArr=("${bamArr2[@]}"); bamArr2=()
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
for ((i=0;i<${#outPfx[@]}; i++)); do
	inBam=${bamArr[$i]}
	outBam="$outDir/${outPfx[$i]}.${step}.tag.bam"
	printMsg "    BAM: $inBam"
	picard SetNmMdAndUqTags -I $inBam -O $outBam --CREATE_INDEX true --CREATE_MD5_FILE true --REFERENCE_SEQUENCE $opt_f -Xmx${opt_m}g -XX:ParallelGCThreads=$opt_t
	if [ $? -ne 0 ] || [ ! -s $outBam ]; then
		bamArr2+=("$outBam")
		printMsg "    BAM: $inBam tagSAM done"
	else
		echo "Error: tagSAM: $inBam"
		exit 100
	fi
done
bamArr=("${bamArr2[@]}"); 



### Step 3: BaseRecalibrator
((step++))
printMsg "Step${step}: VCF"
outDir="$opt_out_dir/${step}.vcf"
vcfOptions=""

if [ ${#vcfArr[@]} -gt 0 ]; then
	for indVcf in ${vcfArr[@]}; do
		vcfOptions=" $vcfOptions --known-sites $indVcf "
	done
else
	if [ -e "bam.list" ]; then
		rm -rf bam.list
	fi
	for ((i=0;i<${#bamArr[@]}; i++)); do
		echo "${bamArr[$i]}" >> bam.list
	done
	
	bcftools mpileup -b bam.list --fasta-ref $opt_f | bcftools call -mv -o merge.raw.vcf
	cat merge.raw.vcf | bgzip > merge.raw.vcf.gz
	tabix -pvcf merge.raw.vcf.gz

#	vcfutils.pl varFilter -D9999 HS0674.var.vcf >HS0674.varFilter.vcf
#		-D INT    maximum read depth [10000000]
#	bgzip HS0674.varFilter.vcf
#	tabix -pvcf HS0674.varFilter.vcf.gz

	vcfOptions=" --known-sites $outDir/merge.raw.vcf "
fi



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
	time $gatk BaseRecalibrator -R $opt_f -I $inBam $vcfOptions -O $outRecal
	
	if [ $? -ne 0 ] || [ ! -s $outRecal ]; then
		printMsg "    BAM: $inBam BaseRecalibrator done"
	else
		echo "Error: BaseRecalibrator: $inBam"
		exit 100
	fi
	recalArr+=("$outRecal")
done
if [ ${#recalArr[@]} -ne ${#bamArr[@]} ] || [ ${#recalArr[@]} -eq 0 ]; then
	echo "Error: BaseRecalibrator error " >&2
	exit 100
fi



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
	time $gatk ApplyBQSR --bqsr-recal-file ${recalArr[$i]} -R $opt_f -I $inBam -O $outbam
	if [ $? -ne 0 ] || [ ! -s $outbam ]; then
		printMsg "    BAM: $inBam ApplyBQSR done"
	else
		echo "Error: ApplyBQSR: $inBam"
		exit 100
	fi
	bamArr2+=("$outbam")
done
bamArr=("${bamArr2[@]}"); 



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
	time gatk --java-options "$opt_jopt" HaplotypeCaller --emit-ref-confidence GVCF -R $opt_f -I $inBam -O $outVcf -bamout $outBam
	if [ $? -ne 0 ] || [ ! -s $outVcf ]; then
		printMsg "    BAM: $inBam HaplotypeCaller done"
	else
		echo "Error: HaplotypeCaller: $inBam"
		exit 100
	fi
	sample_gvcfs=" $sample_gvcfs -V $outVcf "
done



((step++))
printMsg "Step${step}: CombineGVFs"
outDir="$opt_out_dir/${step}.CombineGVFs"
merged_gvcfs="$outDir/merged.gvcf.gz"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
time gatk CombineGVFs -R $opt_f ${sample_gvcfs} -O $merged_gvcfs
if [ $? -ne 0 ] || [ ! -s $merged_gvcfs ]; then
	printMsg "    CombineGVFs done"
else
	echo "Error: CombineGVFs"
	exit 100
fi



((step++))
printMsg "Step${step}: GenotypeGVCFs"
outDir="$opt_out_dir/${step}.GenotypeGVCFs"
merged_vcf="$outDir/merged.vcf.gz"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
time gatk GenotypeGVCFs -R $opt_f -V $merged_vcf -O $merged_vcf
if [ $? -ne 0 ] || [ ! -s $merged_gvcfs ]; then
	printMsg "    GenotypeGVCFs done"
else
	echo "Error: GenotypeGVCFs"
	exit 100
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
gatk VariantRecalibrator -R $opt_f -V $merged_vcf -mode SNP \
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



gatk ApplyVQSR \
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

if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	echo "CMD used: bedGraphToBigWig $opt_bg $opt_fai $opt_o" >&2
	exit 100
fi
