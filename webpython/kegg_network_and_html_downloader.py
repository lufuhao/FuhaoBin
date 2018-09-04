# coding:utf-8
# Email:fanyucai1@126.com
# 2018.3.15
#http://blog.sina.com.cn/s/blog_83f77c940102xig6.html

import subprocess
import os
import sys

if len(sys.argv) !=3:
    sys.stderr.write("\nThis script is used to download KEGG pathway maps.\n")
    sys.stderr.write("usage:python %s hsa hsa/\n\n" %(sys.argv[0]))
    sys.exit(1)

abb=sys.argv[1]
outdir=sys.argv[2]
pathway = "http://rest.kegg.jp/list/pathway/"
pathway+=abb
subprocess.check_call('mkdir -p %s' % (outdir), shell=True)
os.chdir(outdir)

def ko(pathway,abb):
   subprocess.check_call('curl -s %s -o ID.list' % (pathway), shell=True)
   subprocess.check_call('awk \'{print $1}\' ID.list|awk -F":" \'{print $2}\'|sort -u >abb_ID.list && mv abb_ID.list ID.list',shell=True)
   file = open("ID.list", "r")
   for line in file:
       line = line.strip()
       if not os.path.exists('%s.png' %(line)):
           subprocess.check_call('curl http://rest.kegg.jp/get/%s/image -o %s.png' % (line, line),shell=True)
       if not os.path.exists('%s.html' %(line)):
           subprocess.check_call('curl https://www.kegg.jp/kegg-bin/show_pathway?%s -o %s.html' % (line,line), shell=True)
   file.close()

if __name__ =="__main__":
   ko(pathway, abb)
