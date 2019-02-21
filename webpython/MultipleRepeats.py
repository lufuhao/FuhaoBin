#! /usr/bin/env python
import sys, math, os, argparse  
# Citation:
# 10.1534/g3.118.200948 %J G3: Genes|Genomes|Genetics
# Title: Repeats of Unusual Size in Plant Mitochondrial Genomes: Identification, Incidence and Evolution

# Usage: -din directory of files to find repeats in
#        -word word_size  
  
parser = argparse.ArgumentParser(description='Find repeats in a directory of fasta sequence files')  
parser.add_argument('-din', action='store', dest='din', help='Input .fasta directory')  
parser.add_argument('-word', action='store', dest='word', help='Word size for blast')  
results = parser.parse_args()  
din = results.din  
word = results.word  
  
li = os.listdir(din)  
inputs = filter(lambda x: '.fasta' in x, li)  
inputs.sort()  
  
for i in range(len(inputs)):  
    infile = str(inputs[i])  
    os.system("/home/alan/applications/ROUSFinder.py -m "+word+" "+din+infile)
