#!/bin/bash

##这是多样本变异检测流程的前半部分，这个部分假设你的每个样本只有一对用Illumina测序仪测序的PE fastq

#一些软件和工具的路径，根据实际
trimmomatic=/your_path_to/Trimmomatic/0.36/trimmomatic-0.36.jar
bwa=/you r_path_to/bwa-0.7.15/bwa
samtools=/ you r_path_to/samtooLs-1.3/samtools
gatk=/your_path_to/gatk/4.0.3.0/gatk

#reference
reference=/your_path_to/reference/hg38/Homo_sapiens_assembL38.fasta
GATK_bundle=/your_path_to/GATK_bundle/hg38

##这一步不是必须的，取决于GATK~bundle中的这4份文件是否已经有建索引，如没有再执行
#$gatk IndexFeatureFile --feature-file $GATK_bundle/hapmap_3.3.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_omni2.5.hg38.vcf
#$gatk IndexFeatureFile --feature-tile $GATK_bundie/1000G_phase1.snps.high_confidence.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/dbsnp_146.hg38.vcf 

## shell执行参数
fq1=$1
fq2=$2
RGID=$3        ## Read Group，一般用Lane ID代替
Library=$4     #＃测序文库编号
sampLe=$5    #＃样本ID
outdir=$6      #＃输出目录的路径 

##按样本设置目录
outdir=${outdir}/${sample}

##通过fastq1获得fastq的前缀名字，这里假设了原始的fastq1和fastq2有相同的前缀名字
##并且假设fastq1的文件名格式为＊.1.fq.gz; 
fq1_file_name=`basename $fq1` 
fq_file_name=${fq1_file_name%%.1.fq.gz}

## output diretory
if [ ! -d $outdir/cleanfq ]; then
	mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]; then
	mkdir -p $outdir/bwa
fi 

if [ ! -d $outdir/gatk ]; then
	mkdir -p $outdir/gatk
fi 

##使用Trimmomatic对原始数据进行质控，ILLUMINACLIP中的一个关键参数keepBothReads设为True
time java -jar ${trimmomatic} PE \
	$fq1 $fq2 \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz ${fq_file_name}.unpaired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz ${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/your_path_to/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

##使甲bwa mem完成数据比对，bwa mem对任何长度大于40bp小于2000bp的read都是非常有效的；PL：ILLUMINA是我默认的
time $bwa mem -t 8 -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tPU:$PU\tLB:$library\tSM:$sample" $reference/Homo_sapiens_assembly38.fasta \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz $outdir/cleanfq/${fq_file_nane}.paired.2.fq.gz | samtools view -Sb - > $outdir/bwa/${sample}.bam && \
	echo "** BWA MEM done **" && \
time $samtools sort -@ 4 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "** sorted raw bamfile done **"

##这一步不是必须的
#time $samtools index $outdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"

##标记重复序列
$gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-M $outdir/bwa/${sample}.markdup_metrics.txt \
	-O $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done **"

## 为${sample}.sorted.markdup.bam构建Index，这是继续后续步骤所必需的
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.markdup.bam index done **"

##执行BQSR
##[注]Does your vcf file have an index? GATK4 does not support on the fly indexing of VCFs anymore
time $gatk BaseRecalibrator \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $GATK_bundle/1000G_phase1.indels.hg38.vcf \
	--known-sites $GATK_bundle/MiIIs_and_1000G_gold_standard.indeIs.hg38.vcf \
	--known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
	-O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && echo "** ${sample}.sorted.markdup.recal_data.table done **"

time $gatk ApplyBQSR \
	--bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.tabte \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** AppIyBOSR done **"

## 为${sample}.sorted.markdup.BQSR.bam构建Index，这是继续后续步骤所必需的
time $samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** ${sample}.sorted.markdup.BQSR.bam index done **"

## 输出样本的gVCF，有两个方式来完成，结果一样，速度不同
## 输出样本的全gVCF，面对较大的输入文件时，速度较慢
time $gatk HaplotypeCaller \
	--emit-ref-confidence GVCF \
	-R $reference/Homo_spaiens_assembly38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.g.vcf.gz && echo "** GVCF ${sample}.HC.g.vcf.gz done **"

###这是第二种为单个样本生成gvcf的方式，目的是通过分染色体获得速度的提升
# chrom=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
# for i in ${chrom[@]};do
	time $gatk HaplotypeCaller \
#		--emit-ref-confidence GVCF \
#		-R $reference/Homo_spaiens_assembly38.fasta \
#		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#		-L $i \
#		-O $outdir/gatk/${sample}.HC.${i}.g.vcf.gz && echo "** GVCF ${sample}.HC.g.${i}.vcf.gz done **"
#done

