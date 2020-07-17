#!/usr/bin/env python

#https://mp.weixin.qq.com/s/Dhj5P82Xsh_TzJmNxAhYyA
### Requirements
### 	biopython and bcbio-gff
#pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biopython
#pip install -i https://pypi.tuna.tsinghua.edu.cn/simple bcbio-gf
### 参考 https://biopython.org/wiki/GFF_Parsing
### 使用方式
#python convert_gb_to_gff3.py input.gb output.gff



import sys
from Bio import SeqIO
from BCBio import GFF

in_file = sys.argv[1]
out_file = sys.argv[2]

in_handle = open(in_file)
out_handle = open(out_file,'w')

GFF.write(SeqIO.parse(in_handle,'gb'),out_handle)

in_handle.close()
out_handle.close()
