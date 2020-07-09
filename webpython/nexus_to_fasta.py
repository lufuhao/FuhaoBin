#!/usr/bin/env python

#pip install -U dendropy
#
#python nexus_to_fasta.py input_nexus_dir output_fasta_dir
#
#source: https://mp.weixin.qq.com/s/6YwR2n7F-uEo1CD6y9ijZg
#小明的数据分析笔记本


import os
import sys
import dendropy

in_folder = sys.argv[1]
out_folder = sys.argv[2]

os.makedirs(out_folder)

file_num = len(os.listdir(in_folder))
print("Total ", file_num, " nexus files need to be converted")
print("Ready,go!")

for nex_file in os.listdir(in_folder):
    file_path = in_folder + "/" + nex_file
    file_pre = nex_file.split(".")[0]
    nex = dendropy.DnaCharacterMatrix.get(path=file_path,schema='nexus')
    out_file_path = out_folder + "/" + file_pre + ".fasta"
    nex.write(path=out_file_path,schema='fasta')

print("OK")
