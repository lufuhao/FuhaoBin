#!/usr/bin/env python3
#https://blog.csdn.net/qq_26012913/article/details/110390087
# SnpType.py in.snpeff in.gff3 out
#sort -k1 -k2 -V out > snp_type_F.list
#这个脚本是用来看snp哪些是在基因中（分为exon和intron），哪些是在基因的上游、下游
#最终输出结果为：scaffold position 位置 geneid
#输入文件1：snpeff结果文件（/public/home/wangwen_lab/zhangjiexiong/Ageing/SNP/snpeff/snpeffC_F/02.runC/mor_imp_strictestfilter_Snp_AncRef_selected_anno4.vcf）
#输入文件2：gff3文件（/public/home/wangwen_lab/zhangjiexiong/Ageing/SNP/genome.gff3）

import sys,re
file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]
f1 = open(file1,'r')
f2 = open(file2,'r')
f3 = open(file3,'w')
flag = 'fuck'
f1dick = {}
arr1 = []
for i in f1:
        if re.match("#",i):
                continue
        a = i.split("\t")
        if flag != a[0] and flag is not 'fuck':
                f1dick[flag] = arr1
                flag = a[0]
                arr1 = []
        #       print(a[0])
                arr1.append(a[1])
        elif flag is 'fuck':
                flag = a[0]
        else:
                arr1.append(a[1])
f1dick[flag] = arr1
#for i in f1dick:
#        print(i,f1dick[i])

sc = "fuck"
f2dick1 = {}
f2dick2 = {}
f2arr1 = []
for i in f2:
        if re.match("M",i) == None:
                f2dick1[geneid] = f2arr1
                continue
        a = i.split("\t")
        if a[2] == 'gene':
                f2arr1 = []
                f2arr1.append(a[3:5])
                if sc != a[0]:
                        f2dick2[sc] = f2dick1
                        f2dick1 = {}
                        sc = a[0]
                geneid = a[8].split(";")[0].split("=")[1]
        if a[2] == 'exon':
                f2arr1.append(a[3:5])
f2dick1[geneid] = f2arr1
f2dick2[sc] = f2dick1
del f2dick2['fuck']
f2.close()
for i in f1dick.keys():
        for j in f1dick[i]:#取到scaffold内的所有snp位点
                try:
                        taq = 0
                        for k in f2dick2[i].keys():#对gff3文件中所有的scaffold为i的gene进行遍历
                                a = f2dick2[i][k]
                                amax = max(a[0])
                                amin = min(a[0])
                                flag = 0
                                if int(amin) <= int(j) <= int(amax):#判断在gene中后判断在exon还是在intron
                                        for u in a[1:]:
                                                if  int(min(u))<int(j)<int(max(u)):
                                                        flag = "gene_exon"
                                        if flag != "gene_exon":
                                                flag = "gene_intron"
                                elif int(min(a[0]))-2000 <= int(j) <=int(min(a[0])):#判断是不是在gene下游
                                        flag = 'down_stream'
                                elif int(max(a[0]))+2000 >= int(j) >= int(max(a[0])):#判断是不是gene上游
                                        flag = "up_stream"
                                else:
                                        pass
                                if flag != 0:
                                        f3.write(i+"\t"+j+"\t"+k+"\t"+flag+"\t"+max(a[0])+"\t"+min(a[0])+"\n")
                                        taq = 1
                        if taq == 0:
                                f3.write(i+"\t"+j+"\t"+"\t"+"Uninfluncial"+"\n")

                except KeyError:
                        print("KeyError",sys.exc_info())
f3.close()
