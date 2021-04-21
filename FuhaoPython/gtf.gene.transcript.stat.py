#!/usr/bin/python
#https://mp.weixin.qq.com/s/Yd-afFU1FP17BBlqsWknKg

import sys
import re
args = sys.argv

class Genome_info:
      def __init__(self):
          self.chr = ""
          self.start = 0
          self.end = 0

class Gene(Genome_info):
      def __init__(self):
          Genome_info.__init__(self)
          self.orientation = ""
          self.id = ""

class Transcript(Genome_info):
      def __init__(self):
          Genome_info.__init__(self)
          self.id = ""
          self.parent = ""

class Exon(Genome_info):
      def __init__(self):
          Genome_info.__init__(self)
          self.parent = ""

def main(args):
   """
   input gtf file
   """
   list_chr = []
   list_gene = {}
   list_transcript = {}
   list_exon = []
   with open(args[1]) as fp_gtf:
        for line in fp_gtf:
            if line.startswith("#"):
               continue
            lines = line.strip("\n").split("\t")
            chr = lines[0]
            type = lines[2]
            start = int(lines[3])
            end = int(lines[4])
            orientation = lines[6]
            attr = lines[8]
            if not re.search(r'protein_coding',attr):
               continue
            if not chr in list_chr:
               list_chr.append(chr)

            if type == "gene":
               gene = Gene()
               id = re.search(r'gene_id "([^;]+)";?', attr).group(1)
               gene.chr = chr
               gene.start = start
               gene.end = end
               gene.id = id
               gene.orientation = orientation
               list_gene[id] = gene
               print(id)
                           elif type == "transcript":
               transcript = Transcript()
               id = re.search(r'transcript_id "([^;]+)";?',attr).group(1)
               parent = re.search(r'gene_id "([^;]+)";?', attr).group(1)
               if not parent in list_gene:
                  continue
               transcript.chr = chr
               transcript.start = start
               transcript.end = end
               transcript.id = id
               transcript.parent = parent
               list_transcript[id] = transcript
            elif type == "exon":
                 exon = Exon()
                 parent = re.search(r'transcript_id "([^;]+)";?',attr).group(1)
                 if not parent in list_transcript:
                    continue
                 exon.chr = chr
                 exon.start = start
                 exon.end = end
                 exon.parent = parent
                 list_exon.append(exon)

   chr_gene(list_gene)
   gene_len(list_gene)
   gene_transcript(list_transcript)
   transcript_exon(list_exon)
   exon_pos(list_exon)


def chr_gene(list_gene):
    """
    param list_gene:
    """
    print("Distribution of genes on genome")
    count_gene = {}
    for info in list_gene.values():
        chr = info.chr
        if chr in count_gene:
           count_gene[info.chr] += 1
        else:
           count_gene[info.chr] = 1
    with open("chr_gene.txt", 'w') as fp_out:
         for chr,num in count_gene.items():
             print("\t".join([chr,str(num)]) + "\n")
             fp_out.write("\t".join([chr,str(num)]) +"\n")

def gene_len(list_gene):
    """
    param list_gene:
    """
    print("The length of gene")
    with open("gene_len.txt",'w') as fp_out:
         for gene_id,info in list_gene.items():
             len = info.end - info.start +1
             fp_out.write("\t".join([gene_id,str(len)]) + "\n")
             print("\t".join([gene_id,str(len)]) + "\n")


def gene_transcript(list_transcript):
    """
    :param list_transcript:
    """
    print("The distributions of transcripts on genomes")
    count_transcript = {}
    for info in list_transcript.values():
        gene_id = info.parent
        if gene_id in count_transcript:
           count_transcript[gene_id] += 1
        else:
           count_transcript[gene_id] =1
    with open("gene_transcript.txt",'w') as fp_out:
         for gene_id,num in count_transcript.items():
             print("\t".join([gene_id,str(num)]) + "\n")
             fp_out.write("\t".join([gene_id,str(num)]) + "\n")

def transcript_exon(list_exon):
    """
    param list_exon
    """
    print("Summarise of exons within transcript")
    count_exon = {}
    for exon in list_exon:
        transcript_id = exon.parent
        if transcript_id in count_exon:
           count_exon[transcript_id] += 1
        else:
           count_exon[transcript_id] = 1
    with open("transcript_exon.txt",'w') as fp_out:
         for transcript_id, num in count_exon.items():
             print("\t".join([transcript_id,str(num)]) + "\n")
             fp_out.write("\t".join([transcript_id,str(num)]) + "\n")


def exon_pos(list_exon):
    """
    param list_exon
    """
    print("The coordination of exons")
    count_exon = {}
    for exon in list_exon:
        transcript_id = exon.parent
        if transcript_id in count_exon:
           count_exon[transcript_id] += ",%s-%s" % (str(exon.start),str(exon.end))
        else:
           count_exon[transcript_id] = "%s-%s" % (str(exon.start),str(exon.end))
    with open("exon_pos.txt",'w') as fp_out:
         for transcript_id,pos in count_exon.items():
             print("\t".join([transcript_id,pos])+"\n")
             fp_out.write("\t".join([transcript_id,pos]) + "\n")


if __name__ == "__main__":
   main(args)





