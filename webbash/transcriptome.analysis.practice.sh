#!/bin/bash

###转录组差异表达分析小实战
###https://www.cnblogs.com/wangprince2017/p/9937328.html

######### 1.读文献获取数据 ###########
#文献名称：AKAP95 regulates splicing through scaffolding
#RNAs and RNA processing factors

#查找数据：Data availability
#The RIP-seq an RNA-seq data have been deposited in the Gene Expression Omnibus database, with accession code GSE81916. All other data is available from the author upon reasonable request.

获得GSE号：GSE81916

######### 下载测序数据 #########
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81916获取数据信息，并点击网址下方的ftp，下载测序数据
#从https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA323422可知我们需要的mRNA测序编号为SRR3589956到SRR3589962
#通过Apera下载SRR数据，这里以SRR3589956为例：

ascp -T -i /home/anlan/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR358/SRR3589956/SRR3589956.sra ./

######### 转化fastq测序数据 #########
#通过sratoolkit工具将SRR文件转化为fastq格式的测序数据（写了个shell循环）
for i in $(seq 56 62);do
	nohup fastq-dump --split-3  SRR35899${i} &;
done

#通过fastqc对每个fastq文件进行质检，用multiqc查看整体质检报告（对当前目录下的fastq测序结果进行质检，生成每个fq文件的质检报告总multiqc整合后统计查看）

fastqc *.fastq
multiqc ./
#点击这个url可以查看我这个multiqc报告：http://www.bioinfo-scrounger.com/data/multiqc_report.html

#如果有接头或者质量值不达标的需要进行过滤，这次的数据质量都不错，因此直接进行比对即可


######### 序列比对 #########
#安装hisat2软件，下载人类的hiast2索引文件
#hisat2下载并安装：ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip

#下载hisat2的human索引 ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg19.tar.gz
tar zxvf hg19.tar.gz

#用hisat2进行比对，测序数据放在data目录下，索引文件放在reference/index/hisat2/hg19目录下，SRR3589956-SRR3589958为人的测序数据
for i in $(seq 56 58);do 
	hisat2 -p 4 \
		-x ~/reference/index/hisat2/hg19/genome \
		-1 ./data/SRR35899${i}_1.fastq -2 ./data/SRR35899${i}_2.fastq \
		-S SRR35899$i.sam >SRR35899${i}.log;
done

#用samtools将sam文件转化为bam文件，并使用默认排序
for i in $(seq 56 58);do
	samtools sort -@ 5 -o SRR35899${i}.bam SRR35899${i}.sam;
done

######### reads计数 #########
#用htseq对比对产生的bam进行count计数
#htseq安装，使用miniconda，省事！唯一的问题是htseq版本不是最新的，是0.7.2。想要最新版还是要正常安装，可参考http://www.biotrainee.com/thread-1847-1-2.html
conda install -c bioconda htseq

#用htseq将对比后的结果进行计数
for i in $(seq 56 58);do
	htseq-count -f bam -r pos -s no \
	SRR35899${i}.bam ~/reference/genome/hg19/gencode.v26lift37.annotation.gtf \
	1>SRR35899${i}.count 2>SRR35899${i}_htseq.log;
done

#将3个count文件（SRR3589956.count，SRR3589957.count，SRR3589958.count）合并成一个count矩阵，这是就需要脚本来解决这个问题，不然其他方法会稍微麻烦点

#!/usr/bin/perl -w
use strict;

my $path = shift @ARGV;
opendir DIR, $path or die;
my @dir = readdir DIR;

my $header;
my @sample;
my %hash;
foreach my KaTeX parse error: Expected '}', got 'EOF' at end of input: …dir) { if (file =~ /^\w+.*.count/) {
push @sample, $file;
KaTeX parse error: Undefined control sequence: \t at position 12: header .= "\̲t̲file";
open my $fh, fileordie;while(<
fileordie;while(<fh>) {
chomp;
next if ($_ =~ /^\W+/);
my @array = split /\t/, $_;
KaTeX parse error: Expected '}', got 'EOF' at end of input: hash{array[0]} -> {$file} = $array[1];
}
close KaTeX parse error: Expected 'EOF', got '}' at position 9: fh; }̲ } print "header\n";
map{
my $gene = ;print";​print"gene";
foreach my KaTeX parse error: Undefined control sequence: \t at position 33: … print "\̲t̲".hash{KaTeX parse error: Expected 'EOF', got '}' at position 5: gene}̲ -> {file};
}
print “\n”;
}keys %hash;
按照接下来的剧本，应该讲count_matrix文件导入DESeq进行差异表达分析。但是从这篇文章的Bioinformatic analyses部分可以发现，作者的control组的2组数据是来自2个不同的批次（一个是SRR3589956，另外一个来源GSM1095127 in GSE44976），treat组倒是同一个批次（SRR3589957和SRR3589958）。但是对于Mouse cells来说，倒是满足2个control和2个treat都正常来自同个批次，因此打算重新用SRR3589959-SRR3589962重新做个一个count_matrix进行后续差异分析

转录组差异表达分析小实战（二）
Posted: 八月 14, 2017 Under: Transcriptomics By Kai no Comments

差异基因表达分析
我按照前面的流程转录组差异表达分析小实战（一），将小鼠的4个样本又重新跑了一遍，从而获得一个新的count文件：mouse_all_count.txt，有需要的话，可以下载下来进行后续的差异分析。

一般来说，由于普遍认为高通量的read count符合泊松分布，所以一些差异分析的R包都是基于负二项式分布模型的，比如DESeq、DESeq2和edgeR等，所以这3个R包从整体上来说是类似的（但各自标准化算法是不一样的）。

当然还有一个常用的R包则是Limma包，其中的limma-trend和limma-voom都能用来处理RNA-Seq数据（但对应适用的情况不一样）

下面准备适用DESeq2和edgeR两个R包分别对小鼠的count数据进行差异表达分析，然后取两者的结果的交集的基因作为差异表达基因。

DEseq2

library(DESeq2)
##数据预处理
database <- read.table(file = “mouse_all_count.txt”, sep = “\t”, header = T, row.names = 1)
database <- round(as.matrix(database))

##设置分组信息并构建dds对象
condition <- factor(c(rep(“control”,2),rep(“Akap95”,2)), levels = c(“control”, “Akap95”))
coldata <- data.frame(row.names = colnames(database), condition)
dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

##使用DESeq函数估计离散度，然后差异分析获得res对象
dds <- DESeq(dds)
res <- results(dds)

#最后设定阈值，筛选差异基因，导出数据(全部数据。包括标准化后的count数)
res <- res[order(res$padj),]
diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_DESeq2 <- row.names(diff_gene)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by=“row.names”,sort=FALSE)
write.csv(resdata,file = “control_vs_Akap95.csv”,row.names = F)
最终获得572个差异基因（筛选标准为padj < 0.05, |log2FoldChange| > 1）

edgeR

library(edgeR)
##跟DESeq2一样，导入数据，预处理（用了cpm函数）
exprSet<- read.table(file = “mouse_all_count.txt”, sep = “\t”, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
group_list <- factor(c(rep(“control”,2),rep(“Akap95”,2)))
exprSet <- exprSet[rowSums(cpm(exprSet) > 1) >= 2,]

##设置分组信息，并做TMM标准化
exprSet <- DGEList(counts = exprSet, group = group_list)
exprSet <- calcNormFactors(exprSet)

##使用qCML（quantile-adjusted conditional maximum likelihood）估计离散度（只针对单因素实验设计）
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)

##寻找差异gene(这里的exactTest函数还是基于qCML并且只针对单因素实验设计)，然后按照阈值进行筛选即可
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
diff_gene_edgeR <- subset(tTagKaTeX parse error: Expected 'EOF', got '&' at position 19: …le, FDR < 0.05 &̲ (logFC > 1 | l…table,file = “control_vs_Akap95_edgeR.csv”)
最终获得688个差异基因（筛选标准为FDR < 0.05, |log2FC| > 1）

取DESeq2和edgeR两者结果的交集

diff_gene <- diff_gene_DESeq2[diff_gene_DESeq2 %in% diff_gene_edgeR]
最终的差异基因数目为545个

head(diff_gene)
[1] “ENSMUSG00000003309.14” “ENSMUSG00000046323.8” “ENSMUSG00000001123.15”
[4] “ENSMUSG00000023906.2” “ENSMUSG00000044279.15” “ENSMUSG00000018569.12”
其他两个R包（DESeq和limma）就不在这尝试了，我之前做过对于这4个R包的简单使用笔记，可以参考下：
简单使用DESeq做差异分析
简单使用DESeq2/EdgeR做差异分析
简单使用limma做差异分析

GO&&KEGG富集分析
以前一直没有机会用Y叔写的clusterProfiler包，这次正好看说明用一下。

GO富集，加载clusterProfiler包和org.Mm.eg.db包（小鼠嘛），然后将ENSEMBL ID后面的版本号去掉，不然后面不识别这个ID，然后按照clusterProfiler包的教程说明使用函数即可。

library(clusterProfiler)
library(org.Mm.eg.db)

##去除ID的版本号
diff_gene_ENSEMBL <- unlist(lapply(diff_gene, function(x){strsplit(x, “\.”)[[1]][1]}))
##GOid mapping + GO富集
ego <- enrichGO(gene = diff_gene_ENSEMBL, OrgDb = org.Mm.eg.db,
keytype = “ENSEMBL”, ont = “BP”, pAdjustMethod = “BH”,
pvalueCutoff = 0.01, qvalueCutoff = 0.05)
##查看富集结果数据
enrich_go <- as.data.frame(ego)
##作图
barplot(ego, showCategory=10)
dotplot(ego)
enrichMap(ego)
plotGOgraph(ego)
KEGG富集，首先需要将ENSEMBL ID转化为ENTREZ ID，然后使用ENTREZ ID作为kegg id，从而通过enrichKEGG函数从online KEGG上抓取信息，并做富集

library(clusterProfiler)
library(org.Mm.eg.db)

##ID转化
ids <- bitr(diff_gene_ENSEMBL, fromType = “ENSEMBL”, toType = “ENTREZID”, OrgDb = “org.Mm.eg.db”)
kk <- enrichKEGG(gene = ids[,2], organism = “mmu”, keyType = “kegg”,
pvalueCutoff = 0.05, pAdjustMethod = “BH”, qvalueCutoff = 0.1)
##查看富集结果数据
enrich_kegg <- as.data.frame(kk)
##作图
dotplot(kk)
到这里为止，转录组的差异表达分析算是做完了，简单的来说，这个过程就是将reads mapping 到reference上，然后计数获得count数，然后做差异分析，最后来个GO KEGG，over了。。。

对于mapping和计数这两部还有其实还有好多软件，具体可见文献：Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis，有时间可以都尝试下。

至于GO && KEGG这两步，对于人、小鼠等模式物种来说，不考虑方便因素来说，完全可以自己写脚本来完成，数据可以从gene ontology官网下载，然后就是GO id与gene id相互转化了。KEGG 也是一样，也可以用脚本去KEGG网站上自行抓取，先知道URL规则，然后爬数据即可。
