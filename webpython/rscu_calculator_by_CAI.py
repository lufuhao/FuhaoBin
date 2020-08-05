#!/usr/bin/env python

#堆积柱形图（stacked barplot）展示密码子偏向性的RSCU值 
#https://mp.weixin.qq.com/s?__biz=MzI3NzQ3MTcxMg==&mid=2247484705&idx=1&sn=d035611ff4458996dd8cd0204514c3a9&chksm=eb648daedc1304b8887f8a3c9ad2afd403605d1a021df126b07b2cf47b87c717bec1684acca5&scene=126&sessionid=1596416246&key=221452a4d6b5ef37a496b812549c79f2bab325c9c28779bb94bd5009f2721c2d40f32cfa8ed58ec13e399d203c314c21d6d795fc438af217dc28c6a8ac3c1872c2fd47ba449193d4832335335202d106&ascene=1&uin=MTEyOTU0NzU1NA%3D%3D&devicetype=Windows+10+x64&version=62090529&lang=zh_CN&exportkey=AhYlnL7zhxFYsXolvZGkSpw%3D&pass_ticket=BZH8wpKUj3f0pQVRaHxsVom05HgxCY%2FQhlPexfjINdnk1%2Fdm0K4gn3RqM5JVbuN8
#小明的数据分析笔记本

### Step1: calculate RSCU using python CAI
#rscu_calculator_by_CAI.py
### Step2: RSCU stacked barplot using ggplot2
#rscu_ctacked_barploy.by.ggplot2.Rscript


#########requirements################
# CAI: https://github.com/Benjamin-Lee/CodonAdaptationIndex
#     pip install CAI --user
# BioPython
#     pip install biopython --user



######### Load Modules ##############
import sys
from CAI import RSCU
from Bio import SeqIO



######### subfunction ###############
def ArgsParser(argv):
  functionname = 'ArgsParse'
  partdate = '20200803'
   
  try:
    opts, args = getopt.getopt(argv, "hf:o:", ["help", "fasta=", "output="]) 
　　　　　#表示参数选项有：-h, -f, -o, --help, --fasta, --output，它们相互对应；该方法的返回值有两个元素: 第一个是(opt, value)元组的列表，第二个是一般参数列表，包含那些没有 '-' 或 '--' 的参数
  except getopt.GetoptError:
    print('Error: ***.py -f <fasta file> -o <output>')
    print('  or: ***.py --fasta=<fasta file> --output=<partdate>')
    sys.exit(2)
   
  for opt, arg in opts:  #依次获取列表中的元组项
    if opt in ("-h", "--help"):
      print('***.py -f <fasta file> -o <output>')
      print('or: ***.py --fasta=<functionname> --output=<partdate>')
      sys.exit()
    elif opt in ("-i", "--fasta"):
      fasta_file = arg
    elif opt in ("-o", "--output"):
      output = arg
  print('-----------------------------------------------------------------------')
  print(opts) #元组构成的列表
  print(args) #args指的是不用 '-'或 '--'传递的参数，这里没有传递，所以为空
  print('Input Fasta: ', fasta_file)
  print('Ouput', output)
c2aa = {
        'TGT':'Cys',
        'UGU':'Cys',
        'TGC':'Cys',
        'UGC':'Cys',
        'GAT':'Asp',
        'GAU':'Asp',
        'GAC':'Asp',
        'GAC':'Asp',
        'TCT':'Ser',
        'UCU':'Ser',
        'TCG':'Ser',
        'UCG':'Ser',
        'TCA':'Ser',
        'UCA':'Ser',
        'TCC':'Ser',
        'UCC':'Ser',
        'AGC':'Ser',
        'AGC':'Ser',
        'AGT':'Ser',
        'AGU':'Ser',
        'CAA':'Gln',
        'CAA':'Gln',
        'CAG':'Gln',
        'CAG':'Gln',
        'ATG':'Met',
        'AUG':'Met',
        'AAC':'Asn',
        'AAC':'Asn',
        'AAT':'Asn',
        'AAU':'Asn',
        'CCT':'Pro',
        'CCU':'Pro',
        'CCG':'Pro',
        'CCG':'Pro',
        'CCA':'Pro',
        'CCA':'Pro',
        'CCC':'Pro',
        'CCC':'Pro',
        'AAG':'Lys',
        'AAG':'Lys',
        'AAA':'Lys',
        'AAA':'Lys',
        'ACC':'Thr',
        'ACC':'Thr',
        'ACA':'Thr',
        'ACA':'Thr',
        'ACG':'Thr',
        'ACG':'Thr',
        'ACT':'Thr',
        'ACU':'Thr',
        'TTT':'Phe',
        'UUU':'Phe',
        'TTC':'Phe',
        'UUC':'Phe',
        'GCA':'Ala',
        'GCA':'Ala',
        'GCC':'Ala',
        'GCC':'Ala',
        'GCG':'Ala',
        'GCG':'Ala',
        'GCT':'Ala',
        'GCU':'Ala',
        'GGT':'Gly',
        'GGU':'Gly',
        'GGG':'Gly',
        'GGG':'Gly',
        'GGA':'Gly',
        'GGA':'Gly',
        'GGC':'Gly',
        'GGC':'Gly',
        'ATC':'Ile',
        'AUC':'Ile',
        'ATA':'Ile',
        'AUA':'Ile',
        'ATT':'Ile',
        'AUU':'Ile',
        'TTA':'Leu',
        'UUA':'Leu',
        'TTG':'Leu',
        'UUG':'Leu',
        'CTC':'Leu',
        'CUC':'Leu',
        'CTT':'Leu',
        'CUU':'Leu',
        'CTG':'Leu',
        'CUG':'Leu',
        'CTA':'Leu',
        'CUA':'Leu',
        'CAT':'HIS',
        'CAU':'HIS',
        'CAC':'HIS',
        'CAC':'HIS',
        'CGA':'Arg',
        'CGA':'Arg',
        'CGC':'Arg',
        'CGC':'Arg',
        'CGG':'Arg',
        'CGG':'Arg',
        'CGT':'Arg',
        'CGU':'Arg',
        'AGG':'Arg',
        'AGG':'Arg',
        'AGA':'Arg',
        'AGA':'Arg',
        'TGG':'Trp',
        'UGG':'Trp',
        'GTA':'Val',
        'GUA':'Val',
        'GTC':'Val',
        'GUC':'Val',
        'GTG':'Val',
        'GUG':'Val',
        'GTT':'Val',
        'GUU':'Val',
        'GAG':'Glu',
        'GAG':'Glu',
        'GAA':'Glu',
        'GAA':'Glu',
        'TAT':'Tyr',
        'UAU':'Tyr',
        'TAC':'Tyr',
        'UAC':'Tyr',
    }

 



############# MAIN #####################
if __name__ == '__main__':
	ArgsParser(sys.argv[1:])
	seqs = [rec.seq for rec in SeqIO.parse(fasta_file,'fasta')]
rscu = RSCU(seqs)
fw = open(output,'w')
amino_acid = {}
for aa,bb in rscu.items():
    if c2aa[aa] not in amino_acid:
        amino_acid[c2aa[aa]] = 6
    else:
        amino_acid[c2aa[aa]] -= 0.5
    print(aa,c2aa[aa],round(bb,3),amino_acid[c2aa[aa]])
    fw.write(aa+","+c2aa[aa]+","+str(round(bb,3))+","+str(amino_acid[c2aa[aa]])+"\n")

fw.close()
