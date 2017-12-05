#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)
################# help message ######################################
help() {
cat <<HELP

fq2fa fastq_file

Version: 20171205

Descriptions:
    Convert fastq to fasta
    sed -n '1~4s/^@/>/p;2~4p' fq > fa
    
    support .f[ast]q.gz or .f[ast]q

Options:
  -h    Print this help message

Example:
  fq2fa xx.fq > xx.fa

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] && help
#################### Defaults #######################################
fastqfile=$1
echo "$fastqfile" >&2
if [[ "$fastqfile" =~ .+\.[fastqFASTQ]{2,5}\.[gzGZ]{2,2}$ ]]; then
	echo "Info: gzipped fastq" >&2
	zcat $fastqfile | sed -n '1~4s/^@/>/p;2~4p'
elif [[ "$fastqfile" =~ .+\.[fastqFASTQ]{2,5}$ ]]; then
	echo "Info: fastq" >&2
	cat $fastqfile | sed -n '1~4s/^@/>/p;2~4p'
else
	echo "Info: unknown fastq format" >&2
	exit 1
fi
if [ $? -ne 0 ]; then 
	echo "Error: failed" >&2
else
	echo "Info: success" >&2
fi
exit 0
