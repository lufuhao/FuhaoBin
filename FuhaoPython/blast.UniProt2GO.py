#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import gzip
import getopt

################# AUTHORS #########################
#  Fu-Hao Lu
#  Post-Doctoral Scientist in Micheal Bevan laboratory
#  Cell and Developmental Department, John Innes Centre
#  Norwich NR4 7UH, United Kingdom
#  E-mail: Fu-Hao.Lu@jic.ac.uk

def ArgsParser(argv):
	functionname = 'ArgsParse'
	partdate = '20210113'
	USAGE = "\nDescription: Given a BLAST+ outfmt6 file, associate the query IDs to GO terms according to IDmapping file\n\nusage: ***.py -f <fasta file> -m idmapping.tb.gz -o <output>\n***.py --blast=<fasta file> --mapping idmapping.tb.gz --output=<output>\n\nVersion: partdate\n"

	try:
		opts, args = getopt.getopt(argv, "hf:o:", ["help", "blast=", "output="]) 
#表示参数选项有：-h, -f, -o, --help, --fasta, --output，它们相互对应；该方法的返回值有两个元素: 第一个是(opt, value)元组的列表，第二个是一般参数列表，包含那些没有 '-' 或 '--' 的参数
	except getopt.GetoptError:
		print(USAGE)
		print('Error: invalid arguments')
		sys.exit(2)

	for opt, arg in opts:  #依次获取列表中的元组项
		if opt in ("-h", "--help"):
			print(USAGE)
			sys.exit(1)
		elif opt in ("-i", "--blast"):
			fasta_file = arg
		elif opt in ("-o", "--output"):
			output = arg
		elif opt in ("-m", "--mapping"):
			outMapping = arg
	print('-----------------------------------------------------------------------')
	print(opts) #元组构成的列表
	print(args) #args指的是不用 '-'或 '--'传递的参数，这里没有传递，所以为空
	print('Input Fasta: ', fasta_file)
	print('Mapping file: ', outMapping)
	print('Ouput', output)
	return [fasta_file, outMapping, output]



def parseIDmapping(filename):
	UniProt_GO = {}
	with gzip.open(filename, 'r') as f:
		for line in f:
			lsplit = line.rstrip().split("\t")
			if lsplit[7]:
				UniProt_GO[lsplit[1]] = lsplit[7]
	return UniProt_GO

def parseBlastOut(filename):
	tab_res = {}
	with open(filename, 'r') as f:
		for line in f:
			lsplit = line.split()
			tab_res.setdefault(lsplit[0], set()).add(lsplit[1])
	return tab_res


if __name__ == '__main__':
	in1, id1, out2=ArgsParser(sys.argv[1:]) #因为sys.argv[0]是脚本名称
	UniProtKB_GO = parseIDmapping(id1)
	BlastOut = parseBlastOut(in1)
	
	OUT = open(out2, 'w')
	for i in BlastOut:
		temp = []
		for j in BlastOut[i]:
			if j in UniProtKB_GO:
				go = UniProtKB_GO[j].split("; ")
				temp = temp + go
			else:
				continue
		if temp:
			OUT.write(i + "\t" +  ",".join(set(temp)) + "\n")

	OUT.close()

	print('Info: job Done')
	sys.exit(2)
