#!/bin/bash
source FuhaoBash_CmdMod
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
RunPath=$PWD
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
echo "RunDir: $RunPath"

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


opt_g=''
bamArr=()
pfxArr=()
opt_path_gatk='gatk'
opt_out_dir=$PWD
javaOptions="-Xmx4G"
#################### Subfuctions ####################################
for ((i=0;i<${#bamArr[@]};i++)); do
	${bamArr[$i]}=$(readlink -m "${bamArr[$i]}")
done





#################### Command test ###################################
CmdExit 'gatk'



#################### Defaults #######################################
step=0



#################### Input and Output ###############################
opt_g=$(readlink -m "$opt_g")


#################### Main ###########################################


((step++))
printMsg "Step${step}: CreateSequenceDictionary"
if [ ! -s $opt_g.dict ]; then
	gatk CreateSequenceDictionary -R $opt_g
	if [ $? -ne 0 ] || [ ! -s $opt_g.dict ]; then
		echo "Error: gatk CreateSequenceDictionary: $opt_g" >&2
		exit 100
	fi
fi



### GATK 要求read group的格式
### ID = Read group identifier
###　 　每一个read group 独有的ID，每一对reads 均有一个独特的ID，可以自定义命名；
### PL = Platform
### 　　测序平台；ILLUMINA, SOLID, LS454, HELICOS and PACBIO，不区分大小写；
### SM = sample
### 　　reads属于的样品名；SM要设定正确，因为GATK产生的VCF文件也使用这个名字;
### LB = DNA preparation library identifier
### 　　对一个read group的reads进行重复序列标记时，需要使用LB来区分reads来自那条lane;有时候，同一个库可能在不同的lane上完成测序;为了加以区分，
### 　　同一个或不同库只要是在不同的lane产生的reads都要单独给一个ID. 一般无特殊说明，成对儿read属于同一库，可自定义，比如：library1


((step++))
printMsg "Step${step}: MarkDuplicate"
outDir="$opt_out_dir/${step}.MarkDuplicates"
declare -a bamArr2=()
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
for ((i=0;i<${#bamArr[@]};i++)); do
	gatk MarkDuplicates -I ${bamArr[$i]} -O ${pfxArr[$i]}.${step}.dedup.bam -M ${pfxArr[$i]}.${step}.dedup.metrics.txt
	
	samtools index ${pfxArr[$i]}.${step}.dedup.bam
	
	bamArr2+=("$outDir/${pfxArr[$i]}.${step}.dedup.bam")
done
bamArr=("${bamArr2[@]}");





((step++))
printMsg "Step${step}: GVCF"
outDir="$opt_out_dir/${step}.GVCF"
declare -a bamArr2=()
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
inBAMs=""
for ((i=0;i<${#bamArr[@]};i++)); do
	inBAMs= " $inBAMs -I ${bamArr[$i]} "
done
outVcf="${outDir}/merge.${step}.vcf"
gatk --java-options "$javaOptions" HaplotypeCaller $inBAMs -O merge.${step}.gvcf -R $opt_g
gatk GenotypeGVCFs -R $opt_g -V merge.${step}.gvcf -O merge.${step}.vcf



((step++))
printMsg "Step${step}: SNP_InDel"
outDir="$opt_out_dir/${step}.SNP_InDel"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
outSNP="${outDir}/merge.${step}.SNP.vcf"
outInDel="${outDir}/merge.${step}.InDel.vcf"
## 提取SNP
gatk SelectVariants -V $outVcf -O $outSNP --select-type-to-include SNP
## 提取INDEL
gatk SelectVariants -V $outVcf -O $outInDel --select-type-to-include INDEL

##参数
#--select-type-to-include 选择提取的变异类型{NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED}


### 对vcf文件进行过滤
((step++))
printMsg "Step${step}: filter"
outDir="$opt_out_dir/${step}.filter"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
outSNP2fil="${outDir}/merge.${step}.SNP.filter.vcf"
outInDel2fil="${outDir}/merge.${step}.InDel.filter.vcf"
gatk VariantFiltration -O $outSNP2fil -V $outSNP --filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 ||  SOR > 4.0' --filter-name lowQualFilter --cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing
gatk VariantFiltration -O $outInDel2fil -V $outInDel --filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 ||  SOR > 4.0' --filter-name lowQualFilter --cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing
## 参数
#-O 输出filt.vcf文件
#-V 输入vcf文件
#--filter-expression 过滤条件, VCF INFO 信息
#--cluster-window-size 以10个碱基为一个窗口
#--cluster-size 10个碱基为窗口，若存在3以上个则过滤
#--filter-name 被过滤掉的SNP不会删除，而是给一个标签， 比如 Filter
#--missing-values-evaluate-as-failing 当筛选标准比较多的时候，可能有一些位点没有筛选条件当中的一条或几条，例如下面的这个表达式；QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 并不一定所有位点都有这些信息，这种情况下GATK运行的时候会报很多WARNING信息，用这个参数可以把这些缺少某些FLAG的位点也给标记成没有通过筛选的。



##筛选PASS的SNP，INDEL
((step++))
printMsg "Step${step}: pass"
outDir="$opt_out_dir/${step}.pass"
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi
cd $outDir
## 根据FILTER那列信息进行筛选
grep PASS $outSNP2fil >  merge.${step}.SNP.filter.pass.vcf
grep PASS $outInDel2fil >  merge.${step}.InDel.filter.pass.vcf


exit 0
