#!/bin/bash

GenomeFa=$1
OutPfx=$2



if [ ! -s $GenomeFa ]; then
	echo "Error: Genome File not found" >&2
	exit 100
fi
### Chopping sequences. 
if [ ! -s $OutPfx.1.subseq.fa ]; then
	java -jar $NLR_ANNOTATOR_HOME/ChopSequence.jar -i $GenomeFa -o $OutPfx.1.subseq.fa -l 1000000 -p 10000
	if [ $? -ne 0 ] || [ ! -s $OutPfx.1.subseq.fa ]; then
		echo "Error: ChopSequence running error" >&2
		exit 100
	fi
fi
### NLR-Parser
if [ ! -s $OutPfx.2.subseq.fa.nlr.xml ]; then
	java -jar $NLR_ANNOTATOR_HOME/NLR-Parser3.jar -t 1 -y mast -x $NLR_ANNOTATOR_HOME/meme.xml -i $OutPfx.1.subseq.fa -c $OutPfx.2.subseq.fa.nlr.xml
	if [ $? -ne 0 ] || [ ! -s $OutPfx.2.subseq.fa.nlr.xml ]; then
		echo "Error: NLR-Parser.jar running error" >&2
		exit 100
	fi
fi
### NLR-Annotator
if [ ! -s $OutPfx.3.nlr.gff ]; then
	java -jar $NLR_ANNOTATOR_HOME/NLR-Annotator.jar -i $OutPfx.2.subseq.fa.nlr.xml -o $OutPfx.3.nlr.txt -g $OutPfx.3.nlr.gff -b $OutPfx.3.nlr.bed -m $OutPfx.3.nlr.motifs.bed
#-a $OutPfx.3.nlr.nbarkMotifAlignment.fasta -f $OutPfx.3.nlr.genome.fasta -distanceWithinMotifCombination 500 -distanceForElongating 2500 -distanceBetweenMotifCombinations 10000
	if [ $? -ne 0 ] || [ ! -s $OutPfx.3.nlr.gff ]; then
		echo "Error: NLR-Annotator.jar running error" >&2
		exit 100
	fi
fi
