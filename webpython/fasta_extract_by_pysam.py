#!/usr/bin/env python

import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('list',help='seqIDs to be extracted',type=str)
parser.add_argument('fasta',help='input fasta',type=str)

args = parser.parse_args()
list = args.bed
fasta = args.fasta

fa = pysam.FastaFile(fasta)

with open(list,'r') as f:
     for i in f.readlines():
         if i.strip() != '':
            arr = i.strip().split('\t')
            if arr[0] in fa.references:
               seq = fa.fetch(arr[0])
               print(">%s \n %s" % (arr[0],seq))
