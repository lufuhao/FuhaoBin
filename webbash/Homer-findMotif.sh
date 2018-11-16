#!/bin/bash
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' macs_peaks.bed >homer_peaks.bed
# change 
#		chr1	1454086	1454256	MACS_peak_1	59.88 
#to   
#		MACS_peak_1	chr1	1454086	1454256	+
## we don't need ~/biosoft/homer/bin/findMotifsGenome.pl because we have add all of the script to PATH
findMotifsGenome.pl homer_peaks.bed hg19 motifDir -len 8,10,12
annotatePeaks.pl    homer_peaks.bed hg19 1>peakAnn.xls 2>annLog.txt    #-go GODir

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}'  ../macs14_results/highQuaily_peaks.bed >highQuaily_homer.bed
annotatePeaks.pl   highQuaily_homer.bed  hg19 1>highQuaily_peakAnn.xls 2>highQuaily_annLog.txt
awk -F "\t" 'NR>1{ split($8, arr, "("); sub(" ","",arr[1]); print arr[1]}' highQuaily_peakAnn.xls  |sort |uniq -c |awk '{print $2"\t"$1}'
#annotatePeaks.pl peaks.bed hg19 -go GODir 1>peakAnn.xls 2>annLog.txt 

#	make sure that you have install Homer successfully and add the bin folder to the path :
#	vi ~/.bashrc  and add : export PATH="$PATH:/home/jmzeng/biosoft/homer/bin"

# To run Homer, you need to install it first.
# It's a little complicate 

## Download and install blat
cd ~/biosoft
mkdir blat &&  cd blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient   
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer   
cp blat ~/my-bin/bin 

## Download and install Ghostscript
cd ~/biosoft
mkdir Ghostscript &&  cd Ghostscript
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs919/ghostscript-9.19-linux-x86_64.tgz
tar zxvf ghostscript-9.19-linux-x86_64.tgz
cd ghostscript-9.19-linux-x86_64/
cp gs-919-linux_x86_64  ~/my-bin/bin/gs

## Download and install weblogo
cd ~/biosoft
mkdir weblogo &&  cd weblogo
wget http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz
cp seqlogo ~/my-bin/bin/
# WebLogo requires a recent version of gs (ghostscript) 
# <http://www.cs.wisc.edu/~ghost/> to create PNG and PDF output, and 
# ImageMagic's convert <http://www.imagemagick.org/> to create GIFs.
# 
# ./seqlogo -F PDF -f globin.fasta > ../globin.pdf   
# ./seqlogo -F PNG -f globin.fasta > ../globin.png
# ./seqlogo -F GIF -f globin.fasta > ../globin.png
# 

## Download and install seqlogo
# http://bioconductor.org/packages/release/bioc/html/seqLogo.html

## Download and install homer (Hypergeometric Optimization of Motif	EnRichment)
## // http://homer.salk.edu/homer/
## // http://blog.qiubio.com:8080/archives/3024	
## The commands gs, seqlogo, blat, and samtools should now work from the command line
cd ~/biosoft
mkdir homer &&  cd homer
wget http://homer.salk.edu/homer/configureHomer.pl 
perl configureHomer.pl -install
perl configureHomer.pl -install hg19
