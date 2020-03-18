#!/usr/ bin/bash

##这个流程假设你只有一个样本，这个样本只有一对用Illumina测序仪测序的PE fastq数据文件。

##一些软件和工具的路径，根据实际
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


## 对于单个样本来说，有四个完成变异检测的方式，结果是一样的，可以按照需要挑选一种（以下默认第一种）。
## 第一，直接调用HaplotypeCaller输出样本VCF，面对较大的输入文件时，速度较慢
tIme $gatk HaplotypeCaller \
	-R $reference/Homo_sapiens_assembty38.fasta \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"

###第二，输出这个样本每个染色体的vcf，然后再合并所有的染色体结果．目的是提高逗度，这不是必须的．仅笼通过分染色体获得速度的提升
# chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
# for i in ${chrom[@]}; do
#	time $gatk HaplotypeCaller \
#		-R $reference/Homo.sapiens_assembty38.fasta \
#		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#		-L $i \
#		-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done **" &
# done && wait
#merge_vcfs=""
#for i in ${chrom[@]}; do
#	merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sampte}.HC.${i}.vcf.gz \\"\n
#done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

###第三，先输出样本的全gVCF，再进行GenotypeGVCFs，这个方式在单样本情况下不是必须的，但是多样本的标配，面对较大的输入文件时，速度较慢
#time $gatk HaplotypeCaller \
#	--emit-ref-confidence GVCF \
#	-R $reference/Homo_sapiens_assembly38.fasta \
#	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#	-O $outdir/gatk/${sample}.HC.g.vcf.gz && echo "** ${sample}.HC.g.vcf.gz done **"
#time $gatk GenotypeGVCF \
#	-R $reference/Homo.sapiens_assembly38.fasta \
#	-V $outdir/gatk/${sample}.HC.g.vcf.gz \
#	-O $outdir/gatk/${sample).HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done **"
#
###第四，输出每个染色体的gvcf，然后对每个染色体单独进行GenotypeGVCFs，目的是提高速度，这不是必须的，仅是通过分染色体获得速度的提升
# chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
#for i in ${chrom[@]};do
#	time $gatk HaplotypeCaller \
#		--emit-ref-confidence GVCF \
#		-R $reference/Homo_sapiens_assembly38.fasta \
#		-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#		-L $i \
#		-O $outdir/gatk/${sample}.HC.${i}.q.vcf.gz && \
#	time $gatk GenotypeGVCFs \
#		-R $reference/Homo_sapiens_assembly38.fasta \
#		-V $outdir/gatk/${sample}.HC.${i}.g.vcf.gz \
#		-O $outdir/gatk/${sample}.HC.${i}.vcf.qz && echo "** ${sample}.HC.${i}.vcf.gz done **" &
#done && wait
#merge_vcfs=""
#for i in ${chrom[@]};do
#	merge_vcfs=${merge_vcfs}“ -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
#done && time $gatk MerqeVcfs ${merqe_vcfs} -O $outdir/gatk/${sample}.HC.vcf.qz && echo "** MergeVcfs done **"

##VQSR，由于评价SNP和indel质量高低的标准是不同的，因此，需要分SNP和Indel这两种不同的摸式，分别进行质控
##首先是SNP mode
time $gatk VariantRecalibrator \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
	-resource:1000G,known=false,training=true,truth=faise,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/dbsnp_146.hg38.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	-rscriptFile $outdir/gatk/${sample}.HC.snps.plots.R \
	-tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-O $outdir/gatk/${sample}.HC.snps.recal && \
time $gatk ApplyVQSR \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	--ts-filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-recalFile $outdir/gatk/${sample}.HC.snps.recal \
	-mode SNP \
	-O $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"
##然后是Indel mode
time $gatk VariantRecalibrator \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	-resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundte/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	 -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 6 \
	-rscriptFile $outdir/gatk/${sample}.HC.snps.indels.plots.R \
	-tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
	-O $outdir/gatk/${sample}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-input $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	--ts-filter_level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
	-recalFile $outdir/gatk/${sample}.HC.snps.indels.recal \
	-mode INDEL \ 
	-O $outdir/gatk/${sample}.HC.VQSR.vcf.gz && echo "**SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

#可以被删除清理的文件，这不是必须执行的
#1）对于比对文件只有最终的${sample}.sorted.markdup.BQSR.bam值得留下来
#2）对于VCF，可以保留.g.vcf.gz，原始HC.vcf.gz和完成质控的HC.VQSR.vcf.gz
#rm -f $outdir/bwa/${sample}.bam $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.sorted.markdup.bam*
#rm -f $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz && rm -f $outdir/gatk/${sample}.HC.*.recal
