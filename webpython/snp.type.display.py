#!/usr/bin/env python3
#https://blog.csdn.net/qq_26012913/article/details/110390087
#python display_SnpType.py snp_type_F.list snp_type_F_clear.list
#这个是用来处理SnpType.py
import sys,re,math
file1 = sys.argv[1]
file2 = sys.argv[2]
f1 = open(file1,'r')
f2 = open(file2,'w')
la = [0,0]
a = []
fuck  = []
arr =[]
arrline = []
for i in f1:
        if re.search("Uninfluncial",i):
                f2.write(i)
        else:
                line = i.strip()
                a = line.split("\t")
                if a[1] == la[1]:
                        arr.append(la[4])
                        arr.append(la[5])
                        arrline.append(la)
                else:
                        if la[1] == 0:
                                la = a
                                continue
                        arrline.append(la)
                        arr.append(la[4])
                        arr.append(la[5])
                        for j in arr:
                                c = abs(int(j)-int(la[1]))
                                fuck.append(c)
                        t = min(fuck)
                        zjx = fuck.index(t)
                        x = zjx//2
                        fuck = []
                        for zz in arrline[x]:
                                f2.write(zz+"\t")
                        f2.write("\n")
                        arrline = []
                        arr = []
                la = a
arrline.append(la)
arr.append(la[4])
arr.append(la[5])
for j in arr:
        c = abs(int(j)-int(la[1]))
        fuck.append(c)
t = min(fuck)
zjx = fuck.index(t)
x = zjx//2
fuck = []
for zz in arrline[x]:
        f2.write(zz+"\t")
f2.write("\n")
arrline = []
arr = []
f1.close()
f2.close()
