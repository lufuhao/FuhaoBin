#!/bin/bash

# Suppose we have two different assemblies of the same organism, G1.fasta and G2.fasta, and we want to know the which regions in G2 correspond to our interested regions in G1, what should we do? UCSC LiftOver is a great tool to solve this kind of problems

# split G2 into individual chromosome/scaffold (G2.$i.fa) and save them into ./G2_chr direcotry
# I did this step by using a perl script which I wrote, but you can also use UCSC faSplit tool I think.
# Say there are n chromosomes/scaffolds in the G2 assembly, so we will have G2.1.fa, G2.2.fa, ..., G2.n.fa.
 
# Split the G2 chromosomes/scaffolds into 3K chunks and make lift files
mkdir lift
mkdir split
for i in {1..n}
do
  faSplit size ./G2_chr/G2.$i.fa  3000 ./split/G2.$i.split  -lift=./lift/G2.$i.lft -oneFile
done
 
# run blat
mkdir psl
for i in {1..n}
do
  blat G1.fasta ./split/G2.$i.split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=95 -noHead -minScore=100  ./psl/G2.$i.psl
done
 
# Change coordinates of .psl files to parent coordinate system
mkdir liftup
for i in {1..n}
do
    liftUp -pslQ ./liftup/G2.$i.liftup.psl ./lift/G2.$i.lft warn ./psl/G2.$i.psl
done
 
# Make chain files
mkdir chain_raw
for i in {1..n}
do
    axtChain -linearGap=medium  -faQ -faT -psl ./liftup/G2.$i.liftup.psl ./G1.fasta ./G2.fasta ./chain_raw/$i.chain
done
 
# Merge and sort chain files
chainMergeSort ./chain_raw/*.chain | chainSplit chain_split stdin
 
faSize G1.fasta  -detailed >G1.chr_length.txt
faSize G2.fasta  -detailed >G2.chr_length.txt
 
# Make alignment nets from chain files
mkdir net
for i in  ./chain_split/*.chain
do
 tag=${i/\.\/chain_split\//}
 chainNet $i ./G1.chr_length.txt ./G2.chr_length.txt ./net/$tag.net /dev/null
done
 
# Create liftOver chain file
mkdir over
for i in ./chain_split/*.chain
do
  tag=${i/\.\/chain_split\//}
  netChainSubset  ./net/$tag.net $i ./over/$tag.chain
done
 
cat ./over/*.chain >G1_G2.over.chain
 
# Make bed file to report converted coordinates. We can give the coordinates of our query regions (based on G1 assembly) in the input.bed file and liftOver will report the converted coordinates in conversion.bed file.
liftOver input.bed  ./G1_G2.over.chain conversion.bed unMapped
