#!/bin/bash




refFasta=
bamFile=
knownSites=("../Ref/VCF/dbsnp_138.hg19.vcf" "../Ref/VCF/Mills_and_1000G_gold_standard.indels.hg19.vcf" "../Ref/VCF/1000G_phase1.indels.hg19.vcf")
dbSnps="../Ref/VCF/dbsnp_138.hg19.vcf"
#samtools
#gatk


################################################################################
optKnownSites=""
for indSite in ${knownSites[@]}; do
	if [ -s ${indSite} ]; then
		optKnownSites="${optKnownSites} --known-sites ${indSite}"
	else
		echo "Error: invalid known site for gatk" >&2
		exit 100
	fi
done


################################################################################
outPfxTemp=${bamFile##*/}
outPfx=${outPfxTemp%.bam}
bam1Sort=${outPfx}.1.sort.bam
bam2Dedup=${outPfx}.1.sort.dedup.bam
recalTable=${outPfx}.2.recal.table
bam3Recal=${outPfx}.2.recal.bam
bam2Vcf1=${outPfx}.3.snps.indels.raw.vcf
bam2Vcf2=${outPfx}.3.snps.indels.raw.genotype.vcf
bam2Vcf3=${outPfx}.4.genotype.snps.vcf
bam2Vcf3filter=${outPfx}.4.genotype.snps.filter.vcf
bam2Vcf4=${outPfx}.4.genotype.indels.vcf
bam2Vcf4filter=${outPfx}.4.genotype.indels.filter.vcf
bam2vcf5merge=${outPfx}.4.genotype.SNP.InDels.filter.vcf
bam2vcf6pass=${outPfx}.4.genotype.SNP.InDels.filter.pass.vcf



################################################################################
### Step1: 数据准备：
###用Samtools为参考序列创建一个索引，这是为方便GATK能够快速地获取fasta上的任何序列做准备
samtools faidx $refFasta
###生成.dict文件
gatk CreateSequenceDictionary -R $refFasta -O $refFasta.dict
##dic文件的内容：
#@HD	VN:1.5
#@SQ	SN:chr17	LN:81195210	M5:351f64d4f4f9ddd45b35336ad97aa6de	UR:file:/share/disk5/lianm/GATK4_practice/Ref/chr17.fa



### Step2: 排序及标记重复
gatk SortSam -I $bamFile -O $bam1Sort -R $refFasta -SO coordinate --CREATE_INDEX
### --SORT_ORDER,-SO:SortOrder	queryname/coordinate
### --CREATE_INDEX:Boolean 
# 标记重复，这一步会比较费时
gatk MarkDuplicates -I $bam1Sort -O $bam2Dedup -M ${bam2Dedup%.bam}.metrics --CREATE_INDEX



### Step3: 质量值校正
##建立较正模型
gatk BaseRecalibrator -R $refFasta -I $bam2Dedup -O $recalTable ${knownSites[@]}
##参数说明：
#--known-sites:FeatureInput    One or more databases of known polymorphic sites used to exclude regions around known
#                              polymorphisms from analysis.  This argument must be specified at least once. Required.
##质量校正
gatk ApplyBQSR -R $refFasta -I $bam2Dedup -bqsr $recalTable -O $bam3Recal
#参数说明：
#--bqsr-recal-file,-bqsr:File  Input recalibration table for BQSR  Required.



### Step4: SNP、 INDEL位点识别
##生成gvcf文件
gatk HaplotypeCaller -R $refFasta -I $bam3Recal -ERC GVCF --dbsnp $dbSnps -O $bam2Vcf1
# -L chr17:7400000-7800000
##参数说明：
#     --dbsnp,-D:FeatureInput		dbSNP file  Default value: null.
#     
#     --emit-ref-confidence,-ERC:ReferenceConfidenceMode Mode for emitting reference confidence scores  Default value: NONE.
#     				Possible values: {NONE, BP_RESOLUTION, GVCF}
#     
#     --intervals,-L:String		One or more genomic intervals over which to operate  This argument may be specified 0 or more times. Default value: null.
#    -L 规定识别突变位点的区域，如-L chr17:7400000-7800000 只识别17号染色体7400000-7800000 区域的突变位点。 全外显子组分析请用捕获区域bed文件。
#通过gvcf检测变异
gatk GenotypeGVCFs -R $refFasta --dbsnp $dbSnps --variant $bam2Vcf1 -O $bam2Vcf2
#参数说明：
#--variant,-V:String           A VCF file containing variants  Required.



### Step5: 变异位点过滤
#Steps
#        Extract the SNPs from the call set
#        Apply the filter to the SNP call set
#        Extract the Indels from the call set
#        Apply the filter to the Indel call set
#        Combine SNP and indel call set
#        Get passed call set

#提取SNP位点
gatk SelectVariants -R $refFasta -V $bam2Vcf2 --select-type-to-include SNP -O $bam2Vcf3
#    参数说明：
#     --select-type-to-include,-select-type:Type
#                                   Select only a certain type of variants from the input file  This argument may be specified
#                                   0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL,
#                                   SYMBOLIC, MIXED}

#提取INDEL位点
gatk SelectVariants -R $refFasta -V $bam2Vcf2 --select-type-to-include INDEL -O $bam2Vcf4


#SNP位点过滤
gatk VariantFiltration -R $refFasta -V $bam2Vcf3 --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || \
     ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O $bam2Vcf3filter

#    参数说明：
#     --filter-expression,-filter:String
#                                   One or more expression used with INFO fields to filter  This argument may be specified 0
#                                   or more times. Default value: null.
#     --filter-name:String          Names to use for the list of filters  This argument may be specified 0 or more times.
#                                   Default value: null.

##INDEL位点过滤
gatk VariantFiltration -R $refFasta -V $bam2Vcf4 --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < \
     -20.0" --filter-name "INDEL_FILTER" -O $bam2Vcf4filter

##合并过滤后SNP、 INDEL文件
gatk MergeVcfs -I $bam2Vcf3filter -I $bam2Vcf4filter -O $bam2vcf5merge

##提取PASS突变位点
gatk SelectVariants -R $refFasta -V $bam2vcf5merge -O $bam2vcf6pass -select "vc.isNotFiltered()"

#    参数说明：
#     --selectExpressions,-select:String
#                                   One or more criteria to use when selecting the data  This argument may be specified 0 or
#                                   more times. Default value: null.


exit 0
