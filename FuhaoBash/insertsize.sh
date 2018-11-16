#!/bin/sh
#BSUB -a openmpi
#BSUB -J Try6_2
#BSUB -oo %J.20150625.o
#BSUB -e %J.20150625.e
#BSUB -R "rusage[mem=30000]"
#BSUB -q NBI-Prod128
#BSUB -n 10


#20140612
################Configure##############################
###File
#1stlot
#MP40K2R1="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try1_MP40K/1.fastq/853_LIB6484_LDI5272_NoIndex_L001_R1.fastq.gz"
#MP40K2R2="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try1_MP40K/1.fastq/853_LIB6484_LDI5272_NoIndex_L001_R2.fastq.gz"
#2ndlot
#3rdlot_t-b-6
#MP40K2R1="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try3_MP40K/t-b-6/1020_LIB8566_LDI7201_NoIndex_L001_R1.fastq.gz"
#MP40K2R2="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try3_MP40K/t-b-6/1020_LIB8566_LDI7201_NoIndex_L001_R2.fastq.gz"
#3rdlot_t-b-27
#MP40K2R1="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try3_MP40K/t-b-27/1019_LIB8567_LDI7202_NoIndex_L001_R1.fastq.gz"
#MP40K2R2="/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try3_MP40K/t-b-27/1019_LIB8567_LDI7202_NoIndex_L001_R2.fastq.gz"
#4thplot
#MP40K2R1=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try4_MP40K/1094_LIB9374_LDI7704_NoIndex_L001_R1.fastq.gz
#MP40K2R2=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try4_MP40K/1094_LIB9374_LDI7704_NoIndex_L001_R2.fastq.gz
#5th lot
#MP40K5Lib1R1=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try5_MP40K/1.raw_seq/1288_LIB13586_LDI11307_NoIndex_L001_R1.fastq.gz
#MP40K5Lib1R2=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try5_MP40K/1.raw_seq/1288_LIB13586_LDI11307_NoIndex_L001_R2.fastq.gz
#MP40K5Lib2R1=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try5_MP40K/1.raw_seq/1292_LIB13587_LDI11308_NoIndex_L001_R1.fastq.gz
#MP40K5Lib2R2=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try5_MP40K/1.raw_seq/1292_LIB13587_LDI11308_NoIndex_L001_R2.fastq.gz
#MP40K5Lib1R1=/tgac/workarea/collaborators/luf/Try5_MP40K/1.raw_seq/1288_LIB13586_LDI11307_NoIndex_L001_R1.fastq.gz
#MP40K5Lib1R2=/tgac/workarea/collaborators/luf/Try5_MP40K/1.raw_seq/1288_LIB13586_LDI11307_NoIndex_L001_R2.fastq.gz
#MP40K5Lib2R1=/tgac/workarea/collaborators/luf/Try5_MP40K/1.raw_seq/1292_LIB13587_LDI11308_NoIndex_L001_R1.fastq.gz
#MP40K5Lib2R2=/tgac/workarea/collaborators/luf/Try5_MP40K/1.raw_seq/1292_LIB13587_LDI11308_NoIndex_L001_R2.fastq.gz
#MP40Klib1R1=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try6_MK40K/1460_LIB17562_LDI15017_NoIndex_L001_R1.fastq.gz
#MP40Klib1R2=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try6_MK40K/1460_LIB17562_LDI15017_NoIndex_L001_R2.fastq.gz
##ReadLength250
MP40Klib1R1=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try6_MK40K/1465_LIB17561_LDI15016_NoIndex_L001_R1.fastq.gz
MP40Klib1R2=/nbi/group-data/ifs/NBI/Research-Groups/Mike-Beven/Wheat3DL/2.Raw.seq.data/Try6_MK40K/1465_LIB17561_LDI15016_NoIndex_L001_R2.fastq.gz



###Configure#########################################################
RunDir=/usr/users/celldev/luf/temp/wheat/try6_2
ID1='1465_LIB17561_LDI15016'
threads=10
Contaminants=/usr/users/celldev/luf/temp/wheat/hiseq2500.fasta
#REFERENCE=/usr/users/celldev/luf/DB/v443.4.fa
reference=/usr/users/celldev/luf/DB/chr3b.fa
###PATH
FastQC_PATH=$RunDir/0.fastqc
FLASH_PATH=$RunDir/1.flash
CUTADAPT_PATH=$RunDir/2.cutadapt
TRMMOMAT_PATH=$RunDir/3.trimmomatic
FastqJoin_PATH=$RunDir/4.pair
BWA_PATH=$RunDir/5.bwa
SAMTOOLS_PATH=$RunDir/6.samtools
PICARD_HOME=/usr/users/celldev/luf/local/picard/v1.108/x86_64/bin
threads=10
#####################################################################
if [ ! -d $RunDir ]; then
  mkdir -p $RunDir
fi
cd $RunDir



###0. Fastqccd 
if [ ! -d $FastQC_PATH ]; then
  mkdir -p $FastQC_PATH
fi
cd $FastQC_PATH
mkdir -p group nogroup
fastqc $MP40Klib1R1 -o ./group -f fastq -t $threads --noextract
fastqc $MP40Klib1R2 -o ./group -f fastq -t $threads --noextract

fastqc $MP40Klib1R1 -o ./nogroup -f fastq -t $threads --noextract --nogroup
fastqc $MP40Klib1R2 -o ./nogroup -f fastq -t $threads --noextract --nogroup



###1.FLASH
if [ ! -d $FLASH_PATH ]; then
  mkdir -p $FLASH_PATH
fi
cd $FLASH_PATH
flash -m 20 -r 250 -f 800 -s 200 --output-prefix=$ID1 -O -z -t $threads $MP40Klib1R1 $MP40Klib1R2
#[FLASH] ERROR: Read length must be a positive integer!  Please check option -r.
#overlap.extendedFrags.fastq.gzip
#overlap.notCombined_1.fastq.gz
#overlap.notCombined_2.fastq.gz
fastqc $FLASH_PATH/$ID1.extendedFrags.fastq.gz -o . -f fastq -t $threads --noextract --nogroup
fastqc $FLASH_PATH/$ID1.notCombined_1.fastq.gz -o . -f fastq -t $threads --noextract --nogroup
fastqc $FLASH_PATH/$ID1.notCombined_2.fastq.gz -o . -f fastq -t $threads --noextract --nogroup

zcat $ID1.extendedFrags.fastq.gz | perl -ne 'chomp; if ($_=~/^\@(\S+)\s*/) {$ID=$1;}else{die "IDmatcherror\n";} $seq=<>; chomp $seq; $length=length($seq); <>; <>; print $ID."\t".$length."\n";' > $ID1.extendedFrags.fastq.length
perl /usr/users/celldev/luf/lufuhao/SizeCollectBin_luf.pl $ID1.extendedFrags.fastq.length 20 > $ID1.extendedFrags.fasta.SizeBin
cat $ID1.extendedFrags.fastq.length | cut -f 1 > $ID1.extendedFrags.fastq.ID
seqtk subseq $MP40Klib1R1 $ID1.extendedFrags.fastq.ID > $ID1.seqtk.R1.fq
seqtk subseq $MP40Klib1R2 $ID1.extendedFrags.fastq.ID > $ID1.seqtk.R2.fq
fastq_trim_overlap.pl --flashfile $ID1.extendedFrags.fastq.gz --forward $ID1.seqtk.R1.fq --reverse $ID1.seqtk.R2.fq --outR1 $ID1.nomerge.R1.fq --outR2 $ID1.nomerge.R2.fq
fastq2fasta $ID1.extendedFrags.fastq.gz > $ID1.extendedFrags.fasta
cd-hit-est -i $ID1.extendedFrags.fasta -o $ID1.extendedFrags.cdhit -c 0.90 -n 8 -T 0 -r 1 -d 0 -M 30000
deoverlapLib1R1=$FLASH_PATH/$ID1.notCombined_1.fastq.gz
deoverlapLib1R2=$FLASH_PATH/$ID1.notCombined_2.fastq.gz



###2.cutadapt
if [ ! -d $CUTADAPT_PATH ]; then
  mkdir -p $CUTADAPT_PATH
fi
cd $CUTADAPT_PATH
###Try1 
#Adaptor1='GCCAGCGCT'
#Adaptor2='CGCTGGCAG'
###Try2: t-b-3, t-b-27
#1st Adaptor
#Adaptor1='GATCTGCCAGCGCT'
#Adaptor2='AGCGCTGGCAG'
#2nd t-b Adaptor
Adaptor1='GATCTCTACCAGG'
Adaptor2='CCTGGTAGAG'
Kmer01='CTTCTTCTTCTT'
Kmer02='GAAGAAGAAGAA'
Kmer03='CCCCCCCCCC'
Kmer04='CCTCCTCCTCCT'
cutadapt --format=fastq --times=3 --front=$Adaptor1 --front=$Adaptor2 --overlap=8 -o $ID1.cutadapt_1.fq $deoverlapLib1R1
cutadapt --format=fastq --times=3 --front=$Adaptor1 --front=$Adaptor2 --overlap=8 -o $ID1.cutadapt_2.fq $deoverlapLib1R2

fastqc $ID1.cutadapt_1.fq -o . -f fastq -t $threads --noextract --nogroup
fastqc $ID1.cutadapt_2.fq -o . -f fastq -t $threads --noextract --nogroup

CUTADAPT_R1=$CUTADAPT_PATH/$ID1.cutadapt_1.fq
CUTADAPT_R2=$CUTADAPT_PATH/$ID1.cutadapt_2.fq



### 3. trimmomatic
if [ ! -d $TRMMOMAT_PATH ]; then
  mkdir -p $TRMMOMAT_PATH
fi
cd $TRMMOMAT_PATH
java -jar /usr/users/celldev/luf/local/trimmomatic-0.30/trimmomatic-0.30.jar PE -threads $threads -phred33 \
-trimlog $ID1.R1_R2.trimlog \
$CUTADAPT_R1 \
$CUTADAPT_R2 \
$ID1.R1.trim.fq $ID1.R1.unpaired.fq \
$ID1.R2.trim.fq $ID1.R2.unpaired.fq \
ILLUMINACLIP:"$Contaminants":2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:150
mkdir -p group nogroup
fastqc $ID1.R1.trim.fq -o ./nogroup -f fastq -t $threads --noextract --nogroup
fastqc $ID1.R2.trim.fq -o ./nogroup -f fastq -t $threads --noextract --nogroup

fastqc $ID1.R1.trim.fq -o ./group -f fastq -t $threads --noextract
fastqc $ID1.R2.trim.fq -o ./group -f fastq -t $threads --noextract
TRIMMO_R1=$TRMMOMAT_PATH/$ID1.R1.trim.fq
TRIMMO_R2=$TRMMOMAT_PATH/$ID1.R2.trim.fq


### 4. check Pairness
if [ ! -d $FastqJoin_PATH ]; then
  mkdir -p $FastqJoin_PATH
fi
cd $FastqJoin_PATH
fastq_checkid.pl $TRMMOMAT_PATH/$ID1.R1.trim.fq $TRMMOMAT_PATH/$ID1.R2.trim.fq '\@(\S+)\s*\S*'
if [ $? -eq 0 ]; then
  echo "$ID1 paired"
  touch ${ID1}_paired
else
  echo "$ID1 un-paired"
  touch ${ID1}_unpaired
  exit 1
fi



### 5. BWA
if [ ! -d $BWA_PATH ]; then
  mkdir -p $BWA_PATH
fi
cd $BWA_PATH
ln -sf $reference chr3b.fa
bwa index -a is -p ref chr3b.fa
#fastx_reverse_complement -i $PRESEQ_R1 -o $ID.rc_1.fq -Q33
#fastx_reverse_complement -i $PRESEQ_R2 -o $ID.rc_2.fq -Q33
bwa aln -t 30 ref $TRIMMO_R1 > $ID1.1.sai
bwa aln -t 30 ref $TRIMMO_R2 > $ID1.2.sai
bwa sampe -a 65000 ref $ID1.1.sai $ID1.2.sai $TRIMMO_R1 $TRIMMO_R2 | samtools view -bS -F 12 - > $ID1.bam

cd $BWA_PATH
samtools sort $ID1.bam $ID1.sort
pcrdup -counts $ID1.sort.insertsize.count $ID1.sort.bam 
bamutils.pcrdup.insertsize.cluster.pl $ID1.count 20 20 > $ID1.sort.insertsize.rmdup.count
echo "total unique"
tail -n +6 $ID1.sort.insertsize.rmdup.count | wc -l
tail -n +6 $ID1.sort.insertsize.rmdup.count | perl -ne 'chomp;@arr=split(/\t/);print "$arr[2]\t$arr[3]\n";' > $ID1.sort.insertsize.rmdup.count2
tail -n +6 $ID1.sort.insertsize.rmdup.count | cut -f 4 > $ID1.sort.insertsize.rmdup.depth
tail -n +6 $ID1.sort.insertsize.rmdup.count | cut -f 2 > $ID1.sort.insertsize.rmdup.pos
cat $ID1.sort.insertsize.rmdup.depth | perl -ne 'chomp;$sum+=$_;$count++;END {print "AverageDepth: ". $sum/$count ."\n";}'
tail -n +6  $ID1.sort.insertsize.rmdup.count | cut -f 3 > $ID1.sort.insertsize.rmdup.insert

echo "Unique per pos"
sizebin $ID1.sort.insertsize.rmdup.insert 2000
echo "All perl pos"
sizebin $ID1.sort.insertsize.rmdup.count2 2000
echo "duplicates Bin: 2"
sizebin $ID1.sort.insertsize.rmdup.depth 2
echo "duplicates Bin: 5"
sizebin $ID1.sort.insertsize.rmdup.depth 5


3.insertsize.pl -1 $TRMMOMAT_PATH/$ID1.R1.trim.fq -2 $TRMMOMAT_PATH/$ID1.R2.trim.fq -p 6,7 --reference $REFERENCE --minimum_insert 0 --maximum_insert 60000 --min_mapq 5 --picard_memory 20g --fqformat phred33 --threads $threads --prefix $ID1 


cd /usr/users/celldev/luf/temp/wheat/try5/5.bwa/step1_bwa/
bamutils pcrdup -bam 1288_LIB13586_LDI11307_L0_bwaalnQ5.sort.masked.bam -counts 1288_LIB13586_LDI11307_L0_bwaalnQ5.sort.count /usr/users/celldev/luf/temp/wheat/try5/5.bwa/step1_bwa/1288_LIB13586_LDI11307_L0_bwaalnQ5.sort.bam
cd /usr/users/celldev/luf/temp/wheat/try5/5.bwa2/step1_bwa/
bamutils pcrdup -bam 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.masked.bam -counts 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count /usr/users/celldev/luf/temp/wheat/try5/5.bwa2/step1_bwa/1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.bam
bamutils.pcrdup.insertsize.cluster.pl 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count 20 20 > 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count2
cat 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count2 | cut -f 3 > 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count2.insert
sizebin 1292_LIB13587_LDI11308_L0_bwaalnQ5.sort.count2.insert



if [ ! -d $RunDir/6.trimmer ]; then
  mkdir -p $RunDir/6.trimmer
fi
cd $RunDir/6.trimmer
fastx_trimmer -l 200 -i $TRMMOMAT_PATH/$ID1.R1.trim.fq -o $ID1.R1.trim.200.fq
fastx_trimmer -l 200 -i $TRMMOMAT_PATH/$ID1.R2.trim.fq -o $ID1.R2.trim.200.fq

if [ ! -d $RunDir/7.bwa ]; then
  mkdir -p $RunDir/7.bwa
fi
cd $RunDir/7.bwa
perl ~/local/NGSimple/3.insertsize.pl \
-1 $RunDir/6.trimmer/$ID1.R1.trim.200.fq \
-2 $RunDir/6.trimmer/$ID1.R2.trim.200.fq \
-p 6,7 --reference $REFERENCE \
--minimum_insert 0 --maximum_insert 60000 --min_mapq 5 --picard_memory 20g --fqformat phred33 \
--threads $threads \
--prefix $ID1 

#usearch7.0.1090_i86linux32 -id 0.8 -cluster_fast /usr/users/celldev/luf/temp/wheat/try4/flash/overlap.extendedFrags.fa -uc usearch.out
#wcd -N 10 --show_clusters --show_ext -o cluster12.out -d cluster12.dump --word_len 12 --window_len 120 /usr/users/celldev/luf/temp/wheat/try4/flash/overlap.extendedFrags.fa
#fastqc /usr/users/celldev/luf/temp/wheat/try4/flash/overlap.notCombined_1.fastq.gz -o . -f fastq -t $threads --noextract --nogroup
#fastqc /usr/users/celldev/luf/temp/wheat/try4/flash/overlap.notCombined_2.fastq.gz -o . -f fastq -t $threads --noextract --nogroup
#cd-hit-est -i overlap.fa -o overlap.cdhit -c 0.95 -n 8 -T 0 -r 1 -d 0
#gicl.pl --cpus 10 --query overlap.fa -X -mk 
