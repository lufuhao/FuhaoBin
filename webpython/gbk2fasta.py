#!/usr/bin/env python
'''

python gbk2fasta.py

Author:    huagnls
Data:    2015.7.8,
version:
This script was used to get fa from genbank file ;
*.faa = FASTA Amino Acid file
*.ffn = FASTA nucleotide coding regions file
*.fna = FASTA Nucleic Acid file
'''
import sys, os, argparse, os.path
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script was used to get fa from genbank file; *.faa =pep file; *.ffn=cds file; *fna=genome fa file')
parser.add_argument('-i','--gbk',help='Please input genbank  file',required=True)
parser.add_argument('-o','--out_dir',help='Please input  out_put directory path;default cwd',default = os.getcwd(),required=False)
parser.add_argument('-n','--name',required=False,help='Please specify the output file prefix, default input file prefix')
args = parser.parse_args()
dout=''
if os.path.exists(args.out_dir):
    dout=os.path.abspath(args.out_dir)
else:
    os.mkdir(args.out_dir)
    dout=os.path.abspath(args.out_dir)
args.gbk=os.path.abspath(args.gbk)
if args.name == None :
    dir,suffix=os.path.splitext(args.gbk)
    args.name=os.path.basename(dir)


gbk_file = args.gbk
fna_file = dout+'/'+args.name+".ffn"
faa_file = dout+'/'+args.name+".faa"

input = open(gbk_file, "r")
geneNC = open(fna_file, "w")
geneAA = open(faa_file, "w")

sys.stderr.write("%s.faa = pep file\n%s.ffn = CDS  file\n%s.fna = FASTA genome file\n" %(dout+"/"+args.name,dout+"/"+args.name,dout+"/"+args.name))

SeqIO.convert(gbk_file, "genbank", dout+'/'+args.name+".fna", "fasta")

for seq in SeqIO.parse(input, "genbank") :
    print "Dealing with GenBank file of %s, \nOutput: \ngene nc: %s \ngene aa: %s" % (
        seq.id,
        fna_file,
        faa_file)
    for seq_feature in seq.features :
        geneSeq = seq_feature.extract(seq.seq)
        if seq_feature.type == "CDS" :
            assert len(seq_feature.qualifiers['translation']) == 1
            geneAA.write(">%s %s, %s\n%s\n" % (
                seq_feature.qualifiers['locus_tag'][0],
                seq.name,
                seq.description,
                seq_feature.qualifiers['translation'][0]))
            geneNC.write(">%s %s, %s\n%s\n" % (
                seq_feature.qualifiers['locus_tag'][0],
                seq.name,
                seq.description,
                geneSeq ))

input.close()
geneNC.close()
geneAA.close()
